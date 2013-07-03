#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser

import itertools
from scipy import linalg
import pylab as pl
import matplotlib as mpl


import pandas as p
import sklearn.mixture as mixture
from sklearn import preprocessing
import numpy as np


def main(coverage_file,output_base,n_components_range,cv_types,header=None):
    df = p.read_csv(coverage_file,header=header,index_col=0)
    X = df.values
    X = preprocessing.scale(X)

    lowest_bic = np.infty
    bic = []
    convergence = []
    settings = []
    for cv_type in cv_types:
        for n_components in n_components_range:
            sys.stderr.write('Covariance matrix: '+cv_type+', Number of components: '+str(n_components)+'\n')
            # Fit a mixture of gaussians with EM
            gmm = mixture.GMM(n_components=n_components, covariance_type=cv_type)
            gmm.fit(X)
            labels = gmm.predict(X)
            bic.append(gmm.bic(X))
            convergence.append(gmm.converged_)
            class_series = p.Series(labels,index=df.index)
            setting = str(n_components)+"_" + str(cv_type)
            settings.append(setting)
            class_series.to_csv(output_base+"_" +setting+str(gmm.bic(X))+ '.csv')
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm
    
    bic_series = p.Series(bic,index=settings)
    convergence_series = p.Series(convergence,index=settings)
    bic_df = p.DataFrame({'bic': bic_series, 'converged': convergence_series})
    bic_df.to_csv(output_base+"_BIC_results.csv")
    bic = np.array(bic)
    color_iter = itertools.cycle(['k', 'r', 'g', 'b', 'c', 'm', 'y'])
    clf = best_gmm
    bars = []
    # Plot the BIC scores
    spl = pl.subplot(1, 1, 1)
    for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
        xpos = np.array(n_components_range) + .2 * (i - 2)
        bars.append(pl.bar(xpos, bic[i * len(n_components_range):
                                         (i + 1) * len(n_components_range)],
                           width=0.5, color=color))
    pl.xticks(n_components_range)
    pl.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
    pl.title('BIC score per model')
    xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
        .2 * np.floor(bic.argmin() / len(n_components_range))
    pl.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
    spl.set_xlabel('Number of components')
    spl.legend([b[0] for b in bars], cv_types)
    pl.savefig(output_base+"_bic_comparison.png")


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('coverage',
                        help='specify the coverage file')
    parser.add_argument('-o','--output',
                        help='specify the output base file_name, the number of components and cv_type will be added to this file name.')
    parser.add_argument('--n_start', type=int,
                        help='The number of clusters to start at.')
    parser.add_argument('--n_stop', type=int,
                        help='The number of clusters to stop at.')
    parser.add_argument('--n_step', type=int,
                        help='The step size for the number of clusters.')
    parser.add_argument('--full', action='store_true',
                        help='Add this tag to include full covariance matrix')
    parser.add_argument('--diag', action='store_true',
                        help='Add this tag to include diagonal covariance matrix')
    parser.add_argument('--spherical', action='store_true',
                        help='Add this tag to include spherical covariance matrix')
    parser.add_argument('--tied', action='store_true',
                        help='Add this tag to include a tied covariance matrix')


    parser.add_argument('--header', action='store_true',
                        help='Use this tag if header is included in coverage file')
    args = parser.parse_args()

    cv_types = ['full','diag','spherical','tied']

    n_components_range = range(args.n_start,args.n_stop+1,args.n_step)
    used_cv_types=[]
    for cv_type in cv_types:
        if eval('args.'+cv_type):
            used_cv_types.append(cv_type)
    
    if not args.output:
       sys.exit(-1) 
    if args.header:
        header=0
    else:
        header=None

    main(args.coverage,args.output,n_components_range,used_cv_types,header)

