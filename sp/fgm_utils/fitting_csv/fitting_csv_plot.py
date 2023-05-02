from matplotlib import pyplot as plt
import pandas as pd
from typing import List

def parameter_plot(df: pd.DataFrame, title: str):

    fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(6, 8))
    axes[0].plot(df['att_rot'], df[["G1", "G2", "G3"]])
    axes[0].legend(["G1", "G2", "G3"])
    axes[0].set_ylabel("Gain")
    axes[0].set_title(title)

    axes[1].plot(df['att_rot'], df[["th1", "th2", "th3"]])
    axes[1].set_xlim([-1, 358])
    axes[1].legend(["th1", "th2", "th3"])
    axes[1].set_ylabel("theta")

    axes[2].plot(df['att_rot'], df[["ph1", "ph2", "ph3"]])
    axes[2].legend(["ph1", "ph2", "ph3"])
    axes[2].set_ylabel("phi")

    axes[3].plot(df['att_rot'], df[["O1/G1", "O2/G2", "O3/G3"]])
    axes[3].legend(["O1/G1", "O2/G2", "O3/G3"])
    axes[3].set_ylabel("offset")

    axes[4].plot(df['att_rot'], df[["res_med"]])
    axes[4].set_ylabel("res median")

    plt.tight_layout()
    plt.show()
    

def get_csv(path: str) -> pd.DataFrame:

    df = pd.read_csv(path)

    return df

if __name__ == "__main__":
    
    #path = "2021-03-06_15_elb_attloop_Gthphi.csv"
    #path = "2021-03-17_18_elb_attloop_Gthphi.csv"
    path = "2021-03-17_18_elb_attloop2_Gthphi.csv"
    df = get_csv(path)
    #parameter_plot(df, "att_rot", ["G1", "G2", "G3"], "2021-03-06_15-elb")
    parameter_plot(df, "2021-03-17_18_elb")
