### Get raw CH task for each participant and create a new row of data with all the metrics we're most interested in to later merge with my dataframe

#Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import gzip
import base64
import csv
from io import BytesIO
from math import sqrt
import pingouin as pg
import os
from scipy.stats import norm



# Function that that creates long VCH df from a wide df with compressed JSON data in a column "task_data_taskname"
def load_most_recent_csv(path,prefix):
    files = os.listdir(path)

    # get list of full path of files in the directory
    filepaths = [os.path.join(path, file) for file in files if file.endswith('.csv') and prefix in file and "LABELS" not in file]

    # sort the list of files based on last modification time
    filepaths.sort(key=os.path.getmtime,reverse=True)

    # Print the most recent one so we know what we're getting:
    print(f"Loading the most recent CSV file: {filepaths[0]}")

    # First item is the most recent
    dat = pd.read_csv(filepaths[0])
    dat = dat.drop(columns=[x for x in dat.columns if 'Unnamed' in x or 'level_0' in x or 'index' in x])
    return dat


def get_vch_data(dataframe,task_data_taskname):
    participant_dfs = []
    results_dfs = []

    for index, row in dataframe.iterrows():
        # Extract the JSON string from the column
        # try:
            if isinstance(row[task_data_taskname],str):
                record_id = row['record_id']
                compressed_json = row.get(task_data_taskname, "")

                #decompress and decode the JSON
                decoded_data = base64.b64decode(compressed_json)
                decompressed_json = gzip.decompress(decoded_data).decode('utf-8')

                #Parse JSON data
                data = json.loads(decompressed_json)

                #create dataframes of the variables we want for each block then concatenate into one for each participant, add to list of dataframes for later concatenation
                block_dfs = []
                blocks = ['component_1','component_2','component_3','component_4']
                threshold_mapping = {0: 0, 1: 25, 2: 50, 3: 75}
                block_num = 1
                for block in blocks:
                    df1 = pd.DataFrame()
                    df1['vch_response'] = data[block]['response']
                    df1['vch_rt'] = data[block]['responseTime']
                    df1['contrast'] = data[block]['contrasts']
                    df1['vch_block'] = block_num
                    df1['record_id']=record_id
                    df1['rank'] = df1['contrast'].rank(method='dense').astype(int) - 1
                    df1['vch_intensity'] = df1['rank'].map(threshold_mapping)
                    block_dfs.append(df1)
                    block_num+=1
                
                #add all the blocks into one dataframe
                participant_df = pd.concat(block_dfs,ignore_index=True)

                # #replace the contrast in the intensity column with % of threshold
                # participant_df['vch_intensity'] =participant_df['contrast'].replace({data['processedData']['intensities']['c25']:25,
                # data['processedData']['intensities']['c50']:50,
                # data['processedData']['intensities']['c75']:75,
                # data['processedData']['intensities']['c90']:90

                # })

                #add each participant's dataframe into list of dataframes
                participant_dfs.append(participant_df)

                #create dataframe to hold the results we want to merge back
                df_fin = pd.DataFrame(index=[0])

                df_fin['record_id'] = record_id

                #calculate things we want for final analysis
                df_w_response =participant_df[~(participant_df['vch_response'].isna())]
                df_0_trials = df_w_response[df_w_response['vch_intensity']==0]
                df_fin['total_vch_trials_0'] = len(df_0_trials)

                df_fin['total_vch_trials'] = len(df_0_trials[df_0_trials['vch_response']==1])
                df_fin['total_vch_correct_rejects'] = len(df_0_trials[df_0_trials['vch_response']==0])


                df_real_trials = df_w_response[~(df_w_response['vch_intensity']==0)]
                df_fin['vch_hits']=len(df_real_trials[df_real_trials['vch_response']==1])

                # print(df_real_trials[df_real_trials['response'] == 1])
                # print(len(df_real_trials[df_real_trials['response']==1]))
                # print(df_fin['hits'])

                df_fin['vch_misses']=len(df_real_trials[df_real_trials['vch_response']==0])

                df_fin['vch_hit_rate'] = df_fin['vch_hits']/len(df_real_trials)
                df_fin['vch_miss_rate'] = df_fin['vch_misses']/len(df_real_trials)
                df_fin['vch_false_alarm_rate'] = df_fin['total_vch_trials']/df_fin['total_vch_trials_0']
                df_fin["vch_d_prime"] = norm.ppf(df_fin['vch_hit_rate'])-norm.ppf(df_fin['vch_false_alarm_rate'])

                df_fin['total_vch_trials_25'] = len(df_w_response[df_w_response['vch_intensity']==25])
                df_fin['total_vch_trials_50'] = len(df_w_response[df_w_response['vch_intensity']==50])
                df_fin['total_vch_trials_75'] = len(df_w_response[df_w_response['vch_intensity']==75])

                # Finally! Add threshold from data (the intensity for c75)
                df_fin['vch_threshold'] = data['processedData']['intensities']['c75']


                results_dfs.append(df_fin)
        # except Exception as e:
        #     print(f"here's the record giving us {e}: \n\n{row['record_id']}")
    
    # concatenate each participant's full CH data into one ENORMOUS dataframe
    ch_data_all = pd.concat(participant_dfs,ignore_index=True)

    #concatenate each participant's results into one dataframe for merging with the main dataframe
    ch_results = pd.concat(results_dfs,ignore_index=True)  

    return ch_data_all, ch_results


vch_master, vch_results = get_vch_data(df,'task_data_vch_short_psychedelic_bl')

#Add trial 
vch_master['trial'] = vch_master.groupby('record_id').cumcount() + 1

#Add block
vch_master['block'] = ((vch_master['trial'] - 1) // 30).astype(int) + 1

# merge vch_results with the main dataframe for analysis -- keep the other one as vch_master in case we ever want to use it
vch_results['record_id'] = vch_results['record_id'].astype(int)
df = df.merge(vch_results,on='record_id',how='left')



### Add a straightforward threshold interpretation of intensity
# Define the mapping from rank to threshold
threshold_mapping = {0: 0, 1: 25, 2: 50, 3: 75}

# Group by 'record_id' and apply the transformation
vch_master['rank'] = vch_master.groupby('record_id')['vch_intensity'].rank(method='dense').astype(int) - 1

# Now map the ranks to thresholds
vch_master['vch_intensity_threshold'] = vch_master['rank'].map(threshold_mapping)

# Drop the 'rank' column
vch_master.drop(columns=['rank'], inplace=True)

#Add in vch hits for specific intensities
df['vch_hits_25'] = vch_master[vch_master['vch_intensity']==25].groupby('record_id').sum().reset_index()['vch_response']
df['vch_hits_50'] = vch_master[vch_master['vch_intensity']==50].groupby('record_id').sum().reset_index()['vch_response']
df['vch_hits_75'] = vch_master[vch_master['vch_intensity']==75].groupby('record_id').sum().reset_index()['vch_response']



# #Before moving on, drop any infinities that ended up in the dataframe (since one d prime is infinity)
# df.replace([np.inf, -np.inf], np.nan, inplace=True)

# # Drop rows with where d' is now NaN (indicates they had infinity or didn't have ach/vch data)
# df.dropna(subset=['d_prime','vch_d_prime'],inplace=True)

#Calculate percent yes rate for each stimlus intensity
df['vch_bl_yes_0'] = df['total_vch_trials']/df['total_vch_trials_0']
df['vch_bl_yes_25'] = df['vch_hits_25']/df['total_vch_trials_25']
df['vch_bl_yes_50'] = df['vch_hits_50']/df['total_vch_trials_50']
df['vch_bl_yes_75'] = df['vch_hits_75']/df['total_vch_trials_75']