import sys
import logging
import multiprocessing

from sklearn.mixture import GMM

def parallelized_cluster(args):
    c, cv_type,inits,iters,transform_filter= args
    #Run GMM on the pca transform of contigs with kmer count greater
    #than threshold
    gmm = GMM(n_components=c, covariance_type=cv_type, n_init=inits,
              n_iter=iters).fit(transform_filter)
    bic = gmm.bic(transform_filter)
    if gmm.converged_:
        logging.info("Cluster {0} converged".format(c))
    else:
        logging.warning("Cluster {0} did not converge".format(c))
        print >> sys.stderr, "Cluster {0} did not converge".format(c)
    return bic,c, gmm.converged_


def cluster(max_n_processors,cluster_args):
    #This code should be executed by all threads
    if max_n_processors.use_mpi:
        if max_n_processors.rank != 0:
            cluster_args = []
        cluster_args = max_n_processors.comm.bcast(cluster_args, root=0)
        result = map(parallelized_cluster,cluster_args[max_n_processors.rank::max_n_processors.size])
        #Gather all results to root process again
        results = max_n_processors.comm.gather(result, root=0)
        if max_n_processors.rank == 0:
            results = list(chain(*results))
    
    else:
        pool = multiprocessing.Pool(processes=max_n_processors.size)
        results = pool.map(parallelized_cluster,cluster_args)

    if (max_n_processors.use_mpi and max_n_processors.rank==0) or not max_n_processors.use_mpi:
                bics = [(r[0],r[1]) for r in results]
        Output.write_bic(bics)
        min_bic,optimal_c = min(bics,key=lambda x: x[0])
        gmm = GMM(n_components=optimal_c,covariance_type=cv_type,n_init=inits,
                  n_iter=iters).fit(transform_filter)
    
        return results
