import sys
import logging

from sklearn.mixture import GMM, VBGMM, DPGMM

def cluster(args):
    clusterer = args[0]
    stripped_args = args[1:]
    return clusterer.cluster_method(stripped_args)

class Clusterer(object):
    def __init__(self, algorithm="GMM"):
        if algorithm=="GMM":
            self.cluster_method = cluster_gmm
        elif algorithm=="VBGMM":
            self.cluster_method = cluster_vbgmm
        elif algorithm=="DPGMM":
            self.cluster_method = cluster_dpgmm


def cluster_gmm(args):
    c, cv_type, inits, iters, transform_filter, random_seed, alpha= args
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
    

def cluster_vbgmm(args):
    c, cv_type, inits, iters, transform_filter, random_seed, alpha = args
    vbgmm = VBGMM(n_components=c, covariance_type=cv_type, n_iter=iters,
                  random_state=random_seed, alpha=alpha).fit(transform_filter)
    bic = vbgmm.bic(transform_filter)
    if vbgmm.converged_:
        logging.info("Clustering into {0} clusters converged.".format(c))
    else:
        logging.warning(("Clustering into {0} clusters did not "
                         "converge, consider increasing the number "
                         "of iterations.").format(c))
        print >> sys.stderr, "Cluster {0} did not converge".format(c)
    return bic, c, vbgmm.converged_, vbgmm

def cluster_dpgmm(args):
    c, cv_type, inits, iters, transform_filter, random_seed, alpha = args
    dpgmm = DPGMM(n_components=c, covariance_type=cv_type, n_iter=iters,
                  random_state=random_seed, alpha=alpha).fit(transform_filter)
    bic = dpgmm.bic(transform_filter)
    if dpgmm.converged_:
        logging.info("Clustering into {0} clusters converged.".format(c))
    else:
        logging.warning(("Clustering into {0} clusters did not "
                         "converge, consider increasing the number "
                         "of iterations.").format(c))
        print >> sys.stderr, "Cluster {0} did not converge".format(c)
    return bic, c, dpgmm.converged_, dpgmm
    
