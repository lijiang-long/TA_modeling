import pdb

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd

from family_based_model import update_one_generation,calculate_allele_freq,allele_to_genotype_freq

def heatmap():
    ss = []
    ks = []
    cs = []
    iteration = 1000

    for s in np.linspace(0, 1, num=100):
        for k in np.linspace(0, 1, num=100):
            pp,pw,ww = allele_to_genotype_freq(0.5)
            for i in range(iteration):
                pp,pw,ww =update_one_generation(pp,pw,ww,k=k,s=s,h=0.5)
            new_p,_ = calculate_allele_freq(pp,pw,ww)
            ss.append(s)
            ks.append(k)
            cs.append(new_p)

    df = pd.DataFrame(list(zip(ss, ks,cs)),
               columns =['fitness cost', 'outcrossing rate','final allele frequency'])
    # pdb.set_trace()
    df = df.pivot('outcrossing rate', 'fitness cost','final allele frequency')
    df.columns = ["{:.1f}".format(fnum) for fnum in df.columns]
    df.index = ["{:.1f}".format(fnum) for fnum in df.index]
    fig, ax = plt.subplots()
    # im = ax.imshow(df)

    # plt.show()
    sns.heatmap(df,xticklabels=10,yticklabels=10)
    ax.invert_yaxis()


    # ax.xaxis.set_ticks(np.arange(0, 1, 0.2))
    # ax.set_xticks([0,0.1,0.2])
    # ax.yaxis.set_ticks(np.arange(0, 1, 0.2))
    # sns.lineplot(data=df, x="peelzeel frequency", y="allele change", hue="outcrossing rate")
    # plt.plot([0,1],[0,0])
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    fig.savefig('heatmap.png')


if __name__=='__main__':
    heatmap()

