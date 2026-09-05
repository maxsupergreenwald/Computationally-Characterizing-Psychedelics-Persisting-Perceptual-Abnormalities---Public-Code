"""Helper function `get_field_label_dict`."""

from pathlib import Path as _Path

_DATA_DICT_DIR = str(_Path(__file__).resolve().parent.parent)   # modules/

def get_field_label_dict(field_name, download_path=_DATA_DICT_DIR,
                         datadictprefix="redcap_data_dictionary"):
    
    """
    Returns a dictionary mapping field names to their labels from the most recent data dictionary file.
    
    Args:
        field_name (str): The name of the field to look up.
        download_path (str): The path where the data dictionary files are stored.
        datadictprefix (str): The prefix of the data dictionary files.
        
    Returns:
        dict: A dictionary mapping field names to their labels.
    """
    import os
    import pandas as pd
    
    # Find the most recent data dictionary file
    files = [f for f in os.listdir(download_path) if f.startswith(datadictprefix) and f.endswith('.csv')]
    if not files:
        raise FileNotFoundError("No data dictionary files found.")
    
    latest_file = max(files, key=lambda x: os.path.getmtime(os.path.join(download_path, x)))
    
    # Load the data dictionary
    datadict = pd.read_csv(os.path.join(download_path, latest_file))
    
    # Create a dictionary mapping field names to their labels
    redcaplabelstring = datadict[datadict['Variable / Field Name']== field_name]['Choices, Calculations, OR Slider Labels'].values[0]

    if pd.isna(redcaplabelstring):
        print(f"No labels found for field '{field_name}'. Returning empty dictionary.")
        return {}
    else:
        label_dict = {}
        # Split by ' | ' to get each pair
        pairs = redcaplabelstring.split(' | ')
        for pair in pairs:
            # Split only on the first ', ' to separate key and value
            key, value = pair.split(', ', 1)
            label_dict[int(key.strip())] = value.strip()

        return label_dict

__all__ = ["get_field_label_dict"]
