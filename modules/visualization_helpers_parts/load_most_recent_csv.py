"""Helper function `load_most_recent_csv`."""

import os

def load_most_recent_csv(path,prefix):
    files = os.listdir(path)
    # get list of full path of files in the directory
    filepaths = [os.path.join(path, file) for file in files if file.endswith('.csv') and file.startswith(prefix) and "LABELS" not in file]
    # sort the list of files based on last modification time
    filepaths.sort(key=os.path.getmtime,reverse=True)
    # First item is the most recent
    print(f"Most recent file: {filepaths[0]}")
    return filepaths[0]

__all__ = ["load_most_recent_csv"]
