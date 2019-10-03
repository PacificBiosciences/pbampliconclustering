import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

PALETTE='husl'
sns.set_style('whitegrid')

def plotReads(data,clusters):
    nclusters = len(set(clusters))
    hue       = ['Noise' if c==-1 else c for c in clusters]
    ncol = nclusters/10 + 1
    fig,ax = plt.subplots()
    with sns.color_palette(PALETTE, nclusters):
        sns.scatterplot(data.iloc[:,0],data.iloc[:,1],hue=hue,ax=ax)
#        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), ncol=ncol)
        ax.set_xlabel('PCA1')
        ax.set_ylabel('PCA2')
    fig.tight_layout()
    return fig
