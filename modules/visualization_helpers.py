"""Dynamic loader for visualization helper submodules.

This module auto-loads all public helpers from `visualization_helpers_parts/`.
"""

import importlib
import pkgutil
from pathlib import Path

_PARTS_DIR = Path(__file__).with_name("visualization_helpers_parts")
_PACKAGE_NAME = "visualization_helpers_parts"

__all__ = []

for module_info in sorted(pkgutil.iter_modules([str(_PARTS_DIR)]), key=lambda m: m.name):
    module_name = module_info.name
    if module_name.startswith("_"):
        continue

    module = importlib.import_module(f"{_PACKAGE_NAME}.{module_name}")
    public_names = getattr(module, "__all__", [n for n in dir(module) if not n.startswith("_")])

    for name in public_names:
        globals()[name] = getattr(module, name)
        __all__.append(name)

__all__ = sorted(set(__all__))

# ## Final palette override
try:
    import master_config as _vl
    binary_palette = _vl.binary_palette
    __all__.append('binary_palette')
except Exception:
    pass
