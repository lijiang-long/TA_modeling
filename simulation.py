import pdb

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd

from family_based_model import update_one_generation,calculate_allele_freq,allele_to_genotype_freq


sns.set_style("whitegrid")



def trajectory():
    pp,pw,ww = allele_to_genotype_freq(0.5)
    allele_freqs = [0.5]
    for i in range(100):
        pp,pw,ww =update_one_generation(pp,pw,ww,k=1,s=0.3,h=0.5)
        p,q = calculate_allele_freq(pp,pw,ww)
        allele_freqs.append(p)



def allele_frequency_change_plot():
    ks = []
    ps = []
    cs = []
    iteration = 8

    for k in np.linspace(0, 1, num=11):
        for p in np.linspace(0.01, 0.99, num=100):

            pp,pw,ww = allele_to_genotype_freq(p)
            for i in range(iteration):
                pp,pw,ww =update_one_generation(pp,pw,ww,k=k,s=0.35,h=0.5)
            c_p  = calculate_allele_freq(pp,pw,ww)
            pp,pw,ww =update_one_generation(pp,pw,ww,k=k,s=0.35,h=0.5)
            new_p,_ = calculate_allele_freq(pp,pw,ww)
            ks.append(k)
            ps.append(c_p)
            cs.append(new_p - c_p)

    df = pd.DataFrame(list(zip(ks, ps,cs)),
               columns =['outcrossing rate', 'peelzeel frequency','allele change'])

    sns.lineplot(data=df, x="peelzeel frequency", y="allele change", hue="outcrossing rate")
    plt.plot([0,1],[0,0])
    plt.savefig('outcrossing.png')


if __name__=='__main__':
    # allele_frequency_change_plot()
    heatmap()
    # trajectory()
