import pandas as pd

def display_everything(cancel=False):
    if cancel:
        pd.set_option('display.max_rows', 10)
        pd.set_option('display.max_columns', 6)
        pd.set_option('display.width', 80)
        pd.set_option('display.max_colwidth', 50)
    else:
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', None)        
