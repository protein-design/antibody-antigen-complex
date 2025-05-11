'''
draw plots
'''
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pprint import pprint

class PlotAbAg:

    @staticmethod
    def plot_heatmap(sasa, dist, label):
        fig, ax = plt.subplots(2,1, figsize=(10, 10), layout='tight')
        fig.supxlabel('Residue No')
        fig.supylabel(label)
        
        i = 0
        sns.heatmap(sasa, cmap="YlOrRd", vmax=50, ax=ax[i])
        ax[i].set_xticklabels([])
        ax[i].set_xticks([])
        ax[i].set_yticklabels([])
        ax[i].set_yticks([])
        i = 1
        sns.heatmap(dist, cmap="YlOrRd", vmin=-20, vmax=0, ax=ax[i])
        ax[i].set_yticklabels([])
        ax[i].set_yticks([])
        plt.show()