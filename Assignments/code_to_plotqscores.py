import os
import math
import pylab as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches

"""with open("assignment1.fastq") as x:
    lines=x.readlines()
    Qualityscoreline=lines[3::4]
    np.asarray(a)"""
from Bio import SeqIO
def plot_fastq(filename):
        fastq_parser=SeqIO.parse(filename, "fastq")
        res=[]
        c=0
        for record in fastq_parser:
            score=record.letter_annotations["phred_quality"]
            res.append(score)
            c+=1
            if c>10000:
                break
        df = pd.DataFrame(res)
        l = len(df.T)
        fig ,ax=plt.subplots(figsize=(12,5))
        rect = patches.Rectangle((0,0),l,20,linewidth=0,facecolor='r',alpha=.4)
        ax.add_patch(rect)
        rect = patches.Rectangle((0,20),l,8,linewidth=0,facecolor='yellow',alpha=.4)
        ax.add_patch(rect)
        rect = patches.Rectangle((0,28),l,12,linewidth=0,facecolor='g',alpha=.4)
        ax.add_patch(rect)
        df.mean().plot(ax=ax,c='black')
        boxprops = dict(linestyle='-', linewidth=1, color='black')
        df.plot(kind='box', ax=ax, grid=False, showfliers=False, color=dict(boxes='black',whiskers='black')  )
        ax.set_xticks(np.arange(0, l, 5))
        ax.set_xticklabels(np.arange(0, l, 5))
        ax.set_xlabel('position(bp)')
        ax.set_xlim((0,l))
        ax.set_ylim((0,40))
        ax.set_title('per base sequence quality') 
        plt.show()
        return
plot_fastq("assignment1.fastq")
