import sys
import logging

from sklearn.mixture import GMM

def cluster(args):
    c, cv_type,inits,iters,transform_filter,random_seed= args
    #Run GMM on the pca transform of contigs with kmer count greater
    #than threshold
    gmm = GMM(n_components=c, covariance_type=cv_type, n_init=inits,
              n_iter=iters, random_state=random_seed).fit(transform_filter)
    bic = gmm.bic(transform_filter)
    if gmm.converged_:
        logging.info("Clustering into {0} clusters converged.".format(c))
    else:
        logging.warning(("Clustering into {0} clusters did not "
                         "converge, consider increasing the number "
                         "of iterations.").format(c))
        print >> sys.stderr, "Cluster {0} did not converge".format(c)
    return bic, c, gmm.converged_, gmm

