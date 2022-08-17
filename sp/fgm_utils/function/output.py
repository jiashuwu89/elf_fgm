import numpy as np
import pandas as pd

def output_txt(FGM_datetime, B_x, B_y, B_z, title):
    FGM_data= list(zip(FGM_datetime, B_x, B_y, B_z))
    df = pd.DataFrame(FGM_data, columns=['Timestamp','fsp_gei_x','fsp_gei_y','fsp_gei_z'])
    with open(f"fgm_utils/temp/{FGM_datetime[0][0:4]}{FGM_datetime[0][5:7]}{FGM_datetime[0][8:10]}_{FGM_datetime[0][11:13]}{FGM_datetime[0][14:16]}_{FGM_datetime[-1][11:13]}{FGM_datetime[-1][14:16]}_{title}.txt", 'wb') as f:
        np.savetxt(f, df, fmt='%s%15f%15f%15f')
        f.close()
    
    return