import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

sns.set_style('whitegrid')

def plotEPS(data,minReads):
    distances,indices = NearestNeighbors(n_neighbors=minReads)\
                                        .fit(data)\
                                        .kneighbors(data)
    fig,ax = plt.subplots()
    ax.plot(np.sort(distances[:,-1]))
    ax.set_xlabel('Reads')
    ax.set_ylabel('EPS')
    ax.set_title('Optimal EPS = point of max curvature')
    return fig
