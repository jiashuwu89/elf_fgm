import numpy as np
import pandas as pd

def output_txt(FGM_datetime, columns, column_label, title):

    FGM_data= list(zip(FGM_datetime, *columns))
    df = pd.DataFrame(FGM_data, columns=column_label)
    filename = f"fgm_utils/temp/{FGM_datetime[0][0:4]}{FGM_datetime[0][5:7]}{FGM_datetime[0][8:10]}_{FGM_datetime[0][11:13]}{FGM_datetime[0][14:16]}_{FGM_datetime[-1][11:13]}{FGM_datetime[-1][14:16]}_{title}.csv"
    df.to_csv(filename, index=False)
    
    return