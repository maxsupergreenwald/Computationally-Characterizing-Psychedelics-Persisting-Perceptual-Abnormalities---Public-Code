### create_agent.jl
###
### Defines the generative model shared by EVERY analysis in this directory:
### the perceptual/response model for the conditioned-hallucination (CH) task,
### coupled to a premade binary HGF submodel.
###
### Every other script here — the main fit, BMS, prior recovery, posterior
### recovery, and the PPCs — calls create_agent(), so this file is the single
### point at which the model's structure is defined.  A change here changes
### every result in the manuscript.
###
### ── What the model does on each trial ────────────────────────────────────────
###
### The HGF alone learns the probability that a signal is present.  The CH task
### additionally requires a model of how a percept is formed, because what the
### HGF learns from is the participant's own subjective percept rather than a
### veridical observation.  Each trial therefore runs three steps:
###
###   1. LEARN     Feed the participant's PREVIOUS response into the HGF as its
###                observation, updating the learned log-odds of detection.
###                This closed loop — the model learns from its own prior
###                percept — is what distinguishes the CH-task HGF from a
###                standard binary HGF.
###
###   2. PERCEIVE  Blend the top-down HGF prediction with the bottom-up stimulus
###                intensity, weighted by nu (prior_posterior_weight):
###
###                  xbin_pred = sigmoid(xprob_posterior_mean)
###                  belief    = xbin_pred + 1/(1+nu) * (stimulus - xbin_pred)
###                            = nu/(1+nu) * xbin_pred + 1/(1+nu) * stimulus
###
###                  nu -> inf : percept driven purely by the learned prior
###                  nu = 1    : prior and stimulus weighted equally
###                  nu = 0    : percept driven purely by the stimulus
###
###   3. RESPOND   Pass the belief through a unit-square sigmoid whose steepness
###                is set by beta (action_precision):
###
###                  P(yes) = 0.5 + 0.5 * tanh(beta * (belief - 0.5))
###
###                  large beta : near-deterministic given the belief
###                  small beta : responses pulled toward chance (0.5)
###
### ── Architectures ────────────────────────────────────────────────────────────
###
###   create_agent("binary_2level")   xbin <- xprob
###       xprob's volatility is the fixed free parameter omega2.
###       This is the manuscript's primary model.
###
###   create_agent("binary_3level")   xbin <- xprob <- xvol
###       xvol tracks the time-varying volatility of xprob, adding omega3.
###       Retained for model comparison; omega3 proved unrecoverable.
###
### ── Free parameters ──────────────────────────────────────────────────────────
###
###   action_precision        beta    ActionModel-level
###   prior_posterior_weight  nu      ActionModel-level
###   ("xprob", "volatility") omega2  HGF submodel  (symbol :xprob_volatility)
###   ("xvol",  "volatility") omega3  HGF submodel, 3-level only
###                                   (symbol :xvol_volatility)
###
### The values assigned below are DEFAULTS only.  Fitting scripts override them
### with priors; simulation scripts override them with specific values.  They
### matter solely as the starting point for a model that is never sampled.
###
### ── Recorded state ───────────────────────────────────────────────────────────
###
###   belief   the blended perceptual belief from step 2.  Declared as a State so
###            ActionModels retains it in the agent history, which is what makes
###            it available to the trajectory extraction and the PPCs.

################################################################################
# ─── PACKAGES ─────────────────────────────────────────────────────────────────
################################################################################
using ActionModels, HierarchicalGaussianFiltering, LogExpFunctions


################################################################################
# ─── 1. ACTION MODEL: one trial of learn → perceive → respond ──────────────────
# Called once per trial by ActionModels.  Returns the Bernoulli distribution
# over this trial's response, which Turing uses as the likelihood term.
#
# `input` is the stimulus intensity for the current trial (see §4 of
# model_fitting_singleagent_forarray.jl for how that value is derived).
################################################################################
#### Create Model ####
# ActionModels 0.7 API: action model function receives ModelAttributes, not Agent
update_hgf_binary_conditioned_hallucination = function (model_attributes::ActionModels.ModelAttributes, input::T) where {T<:Union{Real,Missing}}

    ## SETUP ##
    #Extract action model parameters from ModelAttributes
    params = load_parameters(model_attributes)
    prior_posterior_weight = params.prior_posterior_weight
    action_precision = params.action_precision

    #Extract the HGF submodel
    hgf = model_attributes.submodel

    # ── Step 1: LEARN from the previous trial ───────────────────────────────
    # The HGF's observation is the participant's own previous response, not the
    # stimulus.  On trial 1 there is no previous action, so the HGF is left at
    # its initial state and only the prediction is used.
    ## UPDATE HGF ##
    #Participant's own previous action (stored via store_action!) is the input to the HGF
    hgf_input = load_actions(model_attributes).action

    #Only update if previous action is not missing
    if !ismissing(hgf_input)
        update_hgf!(hgf, hgf_input)
    end

    # ── Missing-input path ──────────────────────────────────────────────────
    # With no stimulus to blend in, the percept falls back to the pure top-down
    # prediction.  Used when simulating forward past the end of the observed
    # data; not exercised by the manuscript's complete-data fits.
    ## CALCULATE ACTION PROBABILITY ##
    #Handle missing input
    if ismissing(input)
        #Use current prediction when input is missing
        xprob_posterior_mean = get_states(hgf, ("xprob", "posterior_mean"))
        xbin_prediction_mean = 1 / (1 + exp(-xprob_posterior_mean))
        distribution = Distributions.Bernoulli(xbin_prediction_mean)
        update_state!(model_attributes, :belief, xbin_prediction_mean)
        return distribution
    end

    # ── Step 2: PERCEIVE — blend prediction with stimulus, weighted by nu ───
    #The input is the stimulus strength
    stimulus_strength = input

    #Calculate the prediction for the current trial
    xprob_posterior_mean = get_states(hgf, ("xprob", "posterior_mean"))
    xbin_prediction_mean = 1 / (1 + exp(-xprob_posterior_mean)) #NOTE: COUPLING STRENGTH DEACTIVATED HERE

    #Get the belief as a weighting between stimulus strength and the prediction
    belief =
        xbin_prediction_mean +
        1 / (1 + prior_posterior_weight) * (stimulus_strength - xbin_prediction_mean)


    # ── Step 3: RESPOND — unit-square sigmoid governed by beta ──────────────
    # The max() floor and the tanh form both exist to keep the transform inside
    # its domain for parameter values the sampler may propose; tanh is
    # numerically better behaved here than an explicit piecewise sigmoid.
    #Ensure action_precision is positive
    action_precision_safe = max(1e-10, action_precision)

    #Use a safer calculation to avoid domain errors
    action_probability = 0.5 + 0.5 * tanh(action_precision_safe * (belief - 0.5))

    # ── Guard: reject rather than silently clamp ────────────────────────────
    # A probability outside [0,1] means the sampler has reached a parameter
    # region where the model is not defined.  Throwing RejectParameters makes
    # Turing reject that proposal, which is correct inference.  Clamping the
    # value instead would silently alter the likelihood and bias the posterior.
    #If the action probability is not between 0 and 1
    if !(0 <= action_probability <= 1)
        @show xprob_posterior_mean, xbin_prediction_mean, belief, action_probability
        @show action_precision, prior_posterior_weight, hgf.all_nodes["xprob"].parameters.volatility
        #Throw an error that will reject samples when fitted
        throw(
            RejectParameters(
                "With these parameters and inputs, the action probability became $action_probability, which should be between 0 and 1. Try other parameter settings",
            ),
        )
    end

    #Create Bernoulli distribution with action probability
    distribution = Distributions.Bernoulli(action_probability)

    # ── Record the belief so it survives into the agent history ─────────────
    ## FINALIZE ##
    #Save the belief
    update_state!(model_attributes, :belief, belief)

    #Return the action distribution
    return distribution
end


################################################################################
# ─── 2. AGENT CONSTRUCTOR ──────────────────────────────────────────────────────
# Assembles the action model above with a premade binary HGF submodel of the
# requested depth, and returns the ActionModel that every script here fits,
# simulates from, or evaluates.
################################################################################
#### Create ActionModel ####
function create_agent(hgf_string="binary_3level")

    #Set the action model function
    action_model = update_hgf_binary_conditioned_hallucination

    # ── Default parameter values ────────────────────────────────────────────
    # Placeholders, always overridden by the caller's priors or explicit values.
    # They are also the centres of the priors used throughout, so keeping them
    # here documents where those centres come from:
    #   omega2 = -5.1682685,  nu = 0.72646851,  beta = 0.29350739^(-1) ≈ 3.41
    # (beta is stored as an inverse precision, hence the reciprocal.)
    #Set the default parameters
    parameter_defaults = Dict(
        ("xprob", "volatility") => -5.1682685,
        "prior_posterior_weight" => 0.72646851,
        "action_precision" => 0.29350739^(-1),
    )
    if hgf_string == "binary_3level"
        parameter_defaults[("xvol", "volatility")] = -6
    end

    #Set the HGF (HGF parameters are set here; they are accessible via the submodel interface)
    hgf = premade_hgf(hgf_string, parameter_defaults, verbose=false)

    #Create the ActionModel (ActionModels 0.7 API)
    # - action_precision and prior_posterior_weight are ActionModel-level parameters
    # - xprob_volatility (and xvol_volatility for 3level) are submodel (HGF) parameters,
    #   accessible via set_parameters! fallthrough
    # - action = Action(Bernoulli, Float64): Bernoulli for likelihood, Float64 for storage
    #   (data responses may be Float64 0.0/1.0 from CSV)
    am = ActionModel(action_model;
        parameters = (
            action_precision = Parameter(parameter_defaults["action_precision"]),
            prior_posterior_weight = Parameter(parameter_defaults["prior_posterior_weight"]),
        ),
        states = (belief = State(; discrete=false),),  # Float64 state with missing initial value
        actions = (action = Action(Bernoulli, Float64),),
        observations = (stimulus = Observation(),),
        submodel = hgf,
        verbose = false,
    )

    return am
end
