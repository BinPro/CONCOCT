#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import pandas as p

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/nose_tmp_output'
tmp_basename_dir = tmp_dir_path + '/1'
tmp_basename_dir2 = tmp_dir_path + '/2'
tmp_basename_file = tmp_dir_path + '/file'


CWD = os.getcwd()

class TestCMD(object):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()
        os.mkdir(tmp_basename_dir)
        os.chdir(test_dir_path)

    def tearDown(self):
        """remove temporary output files"""
        for d in os.listdir(tmp_dir_path):
            d_path = os.path.join(tmp_dir_path,d)
            try:
                os.remove(d_path)
            except:
                for f in os.listdir(d_path):
                    f_path = os.path.join(d_path,f)
                    os.remove(f_path)
                os.rmdir(d_path)
        assert os.listdir(tmp_dir_path) == []


    def run_command(self,cov_file='coverage',comp_file='composition.fa',
                    tags=[],basename='nose_tmp_output/1'):
        call_string = "CONCOCT test_data/{0} test_data/{1} --basename {2} -c 3,5,1 2> /dev/null".format(cov_file,comp_file,basename)
        for tag in tags:
            call_string += " " + tag
        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                call_string,
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode

    def run_command_mpi(self,cov_file='coverage',comp_file='composition.fa',
                    tags=[],basename='nose_tmp_output/1'):
        call_string = "mpirun -np 8 CONCOCT test_data/{0} test_data/{1} --basename {2} -c 3,5,1 2> /dev/null".format(cov_file,comp_file,basename)
        for tag in tags:
            call_string += " " + tag
        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                call_string,
                shell=True)
            print >> sys.stderr, "You have mpi support"
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode
            print >> sys.stderr, "You do not have mpi support"

    def file_len(self,fh):
        i=0
        with open(fh) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def md5sum(self,fh):
        infile = open("filename", 'rb')
        content = infile.read()
        infile.close()
        m = hashlib.md5() 
        m.update(content)
        return m.hexdigest()
 
    def test_no_errors(self):
        self.run_command()
        assert_equal(self.c,0,
                     msg = "Command exited with nonzero status")
        run_mpi_test = True
        try:
            import mpi4py
        except ImportError:
            run_mpi_test = False
        if run_mpi_test:
            self.run_command_mpi()
            assert_equal(self.c,0,
                         msg = "Command exited with nonzero status")

    def test_directory_creation(self):
        self.run_command()
        assert_true(isdir(tmp_basename_dir),
                    msg = "Temporary directory not created")
        m_time_first = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        
        # Rerun the concoct and see that the directory is overwritten
        self.run_command()
        m_time_second = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        assert_true(m_time_first != m_time_second,
                     msg = "basename dir is not overwritten")
        L = listdir(tmp_dir_path)
        assert_true(len(L) == 1,
                    msg = "Multiple output directories or files was created")

        # File creation
        self.run_command(basename=tmp_basename_file)
        assert_true(isfile(tmp_basename_file+'_clustering.csv'),
                    msg = "Clustering file is not created, when file is used as basename")
        L = listdir(tmp_basename_dir)
        assert_true(len(L) == 16,
                    msg = "Wrong number of output files, observed {0}".format(L))

    def test_output_files_creation(self):
        # dir as basename
        self.run_command()
        d_p = os.path.join(tmp_basename_dir)
        assert_true(
            isfile(d_p+ '/bic.csv'),
            msg='Bic file is not created'
            )
        assert_true(
            isfile(d_p+ '/clustering_gt1000.csv'),
            msg='Large contigs clustering file is not created'
            )
        assert_true(
            isfile(d_p+ '/pca_means_gt1000.csv'),
            msg='Large contigs cluster pca means file is not created'
            )
        assert_true(
            isfile(d_p+ '/pca_variances_gt1000_dim1.csv'),
            msg='Large contigs cluster pca variances file is not created'
            )
        assert_true(
            isfile(d_p+ '/means_gt1000.csv'),
            msg='Large contigs cluster means file is not created'
            )
        
        assert_true(
            isfile(d_p+ '/variance_gt1000_dim1.csv'),
            msg='Large contigs cluster variance file is not created'
            )
        assert_true(
            isfile(d_p+ '/responsibilities.csv'),
            msg='Large contigs responsibilities file is not created'
            )
        assert_true(
            isfile(d_p+ '/clustering.csv'),
            msg='Clustering file is not created'
            )
        assert_true(
            isfile(d_p+ '/PCA_transformed_data_gt1000.csv'),
            msg='PCA file is not created'
            )
        assert_true(
            isfile(d_p+ '/original_data_gt1000.csv'),
            msg='Original data file is not created'
            )
        assert_true(
            isfile(d_p+ '/log.txt'),
            msg='Log file is not created'
            )
        # dir as file
        self.run_command(basename=tmp_basename_file)
        d_p = tmp_basename_file +'_'
        assert_true(
            isfile(d_p +'bic.csv'),
            msg='Bic file is not created'
            )
        assert_true(
            isfile(d_p+ 'clustering_gt1000.csv'),
            msg='Large contigs clustering file is not created'
            )
        assert_true(
            isfile(d_p+ 'pca_means_gt1000.csv'),
            msg='Large contigs cluster means file is not created'
            )
        assert_true(
            isfile(d_p+ 'pca_variances_gt1000_dim1.csv'),
            msg='Large contigs cluster variances file is not created'
            )
        assert_true(
            isfile(d_p+ 'means_gt1000.csv'),
            msg='Large contigs cluster means file is not created'
            )

        assert_true(
            isfile(d_p+ 'variance_gt1000_dim1.csv'),
            msg='Large contigs cluster variance file is not created'
            )
        assert_true(
            isfile(d_p+ 'responsibilities.csv'),
            msg='Large contigs responsibilities file is not created'
            )
        assert_true(
            isfile(d_p+ 'clustering.csv'),
            msg='Clustering file is not created'
            )
        assert_true(
            isfile(d_p+ 'PCA_transformed_data_gt1000.csv'),
            msg='PCA file is not created'
            )
        assert_true(
            isfile(d_p+ 'original_data_gt1000.csv'),
            msg='Original data file is not created'
            )
        assert_true(
            isfile(d_p+ 'log.txt'),
            msg='Log file is not created'
            )
            
    def test_threshold_functionality(self):
        self.run_command()
        d_p = tmp_basename_dir
        od_1 = d_p+'/original_data_gt1000.csv'
        pca_1 = d_p+'/PCA_transformed_data_gt1000.csv'
        var_1 = d_p+'/variance_gt1000_dim1.csv'
        pca_means_1 = d_p+'/pca_means_gt1000.csv'
        pca_variances_1 = d_p+'/pca_variances_gt1000_dim1.csv'
        means_1 = d_p+'/means_gt1000.csv'
        clust_gt_1 = d_p+'/clustering_gt1000.csv'
        clust_1 = d_p+'/clustering.csv'
        odl_1 = self.file_len(od_1)
        varl_1= self.file_len(var_1)
        pca_meansl_1= self.file_len(pca_means_1)
        pca_variancesl_1= self.file_len(pca_variances_1)
        meansl_1= self.file_len(means_1)
        clust_gtl_1= self.file_len(clust_gt_1)
        clustl_1 = self.file_len(clust_1)
        
        pca_df1 = p.io.parsers.read_table(pca_1)
        pca_m1 = pca_df1.to_records()

        self.run_command(comp_file='composition_some_shortened.fa',
                         basename=tmp_basename_dir2+'/')
        d_p2 = tmp_basename_dir2
        od_2 = d_p2+'/original_data_gt1000.csv'
        pca_2 = d_p2+'/PCA_transformed_data_gt1000.csv'
        var_2 = d_p2+'/variance_gt1000_dim1.csv'
        pca_means_2 = d_p2+'/pca_means_gt1000.csv'
        pca_variances_2 = d_p2+'/pca_variances_gt1000_dim1.csv'
        means_2 = d_p2+'/means_gt1000.csv'
        clust_gt_2 = d_p2+'/clustering_gt1000.csv'
        clust_2 = d_p2+'/clustering.csv'
        odl_2 = self.file_len(od_2)
        varl_2= self.file_len(var_2)
        pca_meansl_2= self.file_len(pca_means_2)
        pca_variancesl_2= self.file_len(pca_variances_2)
        meansl_2= self.file_len(means_2)
        clust_gtl_2= self.file_len(clust_gt_2)
        clustl_2 = self.file_len(clust_2)
        pca_df2 = p.io.parsers.read_table(pca_2)
        pca_m2 = pca_df2.to_records()
        
        assert_true(odl_1!=odl_2,
                    msg='Original data have the same lengths')
        assert_true(varl_1==varl_2,
                    msg='Variance files do not have the same lengths')
        assert_true(pca_meansl_1==pca_meansl_2,
                    msg='PCA mean files do not have the same lengths')
        assert_true(pca_variancesl_1==pca_variancesl_2,
                    msg='PCA variances files do not have the same lengths')
        assert_true(meansl_1==meansl_2,
                    msg='Means files do not have the same lengths')
        assert_true(clust_gtl_1!=clust_gtl_2,
                    msg='Filtered clustering files have the same lengths')
        assert_true(pca_m1.shape!=pca_m2.shape,
                    msg='PCA transformed data has the same shapes')
        assert_true(clustl_1==clustl_2,
                    msg='Clustering files does not have the same lengths')
        assert_true(clust_gtl_2!=clustl_2,
                    msg='Filtered clustering file and full have the same lengths')

    def test_piping_functionality(self):
        f1 = tmp_basename_dir+'/stdout_capture'
        self.run_command(tags=['--pipe','> {0}/stdout_capture'.format(tmp_basename_dir)])

        f1 = open(tmp_basename_dir+'/stdout_capture','rb')
        f1_content = f1.read()
        f1.close()
        f2 = open(tmp_basename_dir + '/clustering.csv','rb')
        f2_content = f2.read()
        f2.close()
        assert_true(len(f1_content)==len(f2_content),
                    msg='stdout and clustering file is not equal')

    def test_bic_sorted(self):
        self.run_command()
        bic = p.io.parsers.read_table(tmp_basename_dir+'/bic.csv',sep=',',index_col=0,header=None)
        assert_true(max(bic.index) == 5,
                    msg='BIC columns are probably mixed up')
        index_l = list(bic.index)
        assert_true(index_l==[3,4,5],
                    msg='BIC file is not sorted')
        # Run command again, to see that file is overwritten
        self.run_command()
        bic = p.io.parsers.read_table(tmp_basename_dir+'/bic.csv',sep=',',index_col=0,header=None)
        assert_true(len(bic)==3,
                    'BIC file is probably appended to')

    def test_logging(self):
        self.run_command()
        with open(tmp_basename_dir+'/log.txt','r') as log:
            log_content = log.read()
            assert_true(len(log_content)>10,
                        "Log content is too small")

    def test_seed(self):
        #Test default behaviour, seed = 11
        first_file = None
        second_file= None
        first_time = None
        second_time = None        

        #Should both run with seed 11
        self.run_command()
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        with open(tmp_basename_dir+'/clustering.csv','r') as clustering:
            first_file=clustering.read()
      
        self.run_command()
        second_time = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        with open(tmp_basename_dir+'/clustering.csv','r') as clustering:
            second_file=clustering.read()
        assert_true(not (first_time==second_time),
                    msg='clustering.csv did not change')
#        assert_true(first_file == second_file,
#                    msg='Clustering outcomes were not the same with same seeds')

        #Should be equal to both above since default seed is 11
	self.run_command(tags=["-f","11"])
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        with open(tmp_basename_dir+'/clustering.csv','r') as clustering:
            first_file=clustering.read()        
        assert_true(not (first_time==second_time),
                    msg='clustering.csv did not change')
        assert_true(first_file == second_file,
                    msg='Clustering outcomes were not the same with same seeds')

        #Test that 0 gives random seed
        first_file = None
        second_file= None
        first_time = None
        second_time = None
        
        #Should give random clustering
        self.run_command(tags=['-f','0'])
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        with open(tmp_basename_dir+'/clustering.csv','r') as clustering:
            first_file=clustering.read()
        
        #Should give random clustering
        self.run_command(tags=['-f','0'])
        second_time = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        with open(tmp_basename_dir+'/clustering.csv','r') as clustering:
            second_file=clustering.read()
        assert_true(not (first_time==second_time),
                    msg='clustering.csv did not change')
        assert_true(first_file == second_file,
                    msg='Clustering outcomes were the same with the different seeds')


        #Test that two differnet seeds give different clustering
        first_file = None
        second_file= None
        first_time = None
        second_time = None
        
        #Should give clustering 2
        self.run_command(tags=['-f','2'])
        first_time = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        with open(tmp_basename_dir+'/clustering.csv','r') as clustering:
            first_file=clustering.read()
        
        #Should give clustering 3
        self.run_command(tags=['-f','3'])
        second_time = os.path.getmtime(tmp_basename_dir+'/clustering.csv')
        with open(tmp_basename_dir+'/clustering.csv','r') as clustering:
            second_file=clustering.read()
        assert_true(not (first_time==second_time),
                    msg='clustering.csv did not change')
        assert_true(not (first_file == second_file),
                    msg='Clustering outcomes were the same with the different seeds')

    def test_log_coverage(self):
        self.run_command()
        original_coverage_data_path = os.path.join(tmp_basename_dir,'original_data_gt1000.csv')
        df = p.io.parsers.read_table(original_coverage_data_path,index_col=0,sep=',')
        # true pseudo count found by running once and using that. Previous tests
        # when we did not normalize the coverage were calculated by hand and 
        # were correct so I assume this is correct with normalization
        true_pseudo_cov = -1.3062 
        calc_pseudo_cov = df.sample_1[0]
        assert_almost_equal(true_pseudo_cov,calc_pseudo_cov,places=4)

    def test_log_coverage_no_cov_normalization(self):
        self.run_command(tags=["--no_cov_normalization"])
        original_coverage_data_path = os.path.join(tmp_basename_dir,'original_data_gt1000.csv')
        df = p.io.parsers.read_table(original_coverage_data_path,index_col=0,sep=',')
        # Manually calculated pseudo coverage using
        # coverage 0.153531, contig length 10132
        true_pseudo_cov = -1.8115 
        calc_pseudo_cov = df.sample_1[0]
        assert_almost_equal(true_pseudo_cov,calc_pseudo_cov,places=4)

    def test_split_pca(self):
        self.run_command(tags=['--split_pca',
                               '--composition_percentage_pca 90',
                               '--coverage_percentage_pca 60'])
        assert_equal(self.c,0,
                     msg = "Split PCA tag results in error")

        d_p = os.path.join(tmp_basename_dir)
        assert_true(
            isfile(d_p+ '/bic.csv'),
            msg='Bic file is not created'
            )
        assert_true(
            isfile(d_p+ '/clustering_gt1000.csv'),
            msg='Large contigs clustering file is not created'
            )
        assert_true(
            isfile(d_p+ '/pca_means_gt1000.csv'),
            msg='Large contigs cluster pca means file is not created'
            )
        assert_true(
            isfile(d_p+ '/pca_variances_gt1000_dim1.csv'),
            msg='Large contigs cluster pca variances file is not created'
            )
        assert_true(
            isfile(d_p+ '/responsibilities.csv'),
            msg='Large contigs responsibilities file is not created'
            )
        assert_true(
            isfile(d_p+ '/clustering.csv'),
            msg='Clustering file is not created'
            )
        assert_true(
            isfile(d_p+ '/PCA_transformed_data_gt1000.csv'),
            msg='PCA file is not created'
            )
        assert_true(
            isfile(d_p+ '/original_data_gt1000.csv'),
            msg='Original data file is not created'
            )
        assert_true(
            isfile(d_p+ '/log.txt'),
            msg='Log file is not created'
            )
        self.run_command(basename=tmp_basename_dir2+'/')
        with open(tmp_basename_dir+'/pca_means_gt1000.csv','r') as pca_m1:
            pca_means1 = pca_m1.read()
        with open(tmp_basename_dir2+'/pca_means_gt1000.csv','r') as pca_m2:
            pca_means2 = pca_m2.read()
        assert_true(not (pca_means1 == pca_means2),
                    msg=('Pca mean files same even with different '
                         'percentage explained variance'))

    def test_versus_reference_results(self):
        self.run_command(tags=["-i 100","-m 1"])
        reference_result_dir = os.path.join(test_dir_path,
                                            'test_data',
                                            'reference_result')
        for ref_f in listdir(reference_result_dir):
            fn = os.path.basename(ref_f)
            # Log will have time stamps and thus not the same
            if fn == 'log.txt':
                continue
            new_f = os.path.join(tmp_basename_dir,fn)
            ref_f = os.path.join(reference_result_dir,fn)
            with open(ref_f,'r') as ref_fh:
                ref_f = ref_fh.read()
            with open(new_f,'r') as new_fh:
                new_f = new_fh.read()
            assert_true(ref_f == new_f,
                        msg=('File not consistent with '
                             'reference file {0}').format(fn))

    def test_versus_reference_results_split_pca(self):
        self.run_command(tags=["--split_pca", "-i 100", "-m 1"])
        reference_result_dir = os.path.join(test_dir_path,
                                            'test_data',
                                            'reference_result_split_pca')
        for ref_f in listdir(reference_result_dir):
            fn = os.path.basename(ref_f)
            # Log will have time stamps and thus not the same
            if fn == 'log.txt':
                continue
            new_f = os.path.join(tmp_basename_dir,fn)
            ref_f = os.path.join(reference_result_dir,fn)
            with open(ref_f,'r') as ref_fh:
                ref_f = ref_fh.read()
            with open(new_f,'r') as new_fh:
                new_f = new_fh.read()
            assert_true(ref_f == new_f,
                        msg=('File not consistent with '
                             'reference file {0}').format(fn))

    def test_versus_reference_results_no_cov_normalization(self):
        self.run_command(tags=["-i 100","-m 1","--no_cov_normalization"])
        reference_result_dir = os.path.join(test_dir_path,
                                            'test_data',
                                            'reference_result_no_cov_normalization')
        for ref_f in listdir(reference_result_dir):
            fn = os.path.basename(ref_f)
            # Log will have time stamps and thus not the same
            if fn == 'log.txt':
                continue
            new_f = os.path.join(tmp_basename_dir,fn)
            ref_f = os.path.join(reference_result_dir,fn)
            with open(ref_f,'r') as ref_fh:
                ref_f = ref_fh.read()
            with open(new_f,'r') as new_fh:
                new_f = new_fh.read()
            assert_true(ref_f == new_f,
                        msg=('File not consistent with '
                             'reference file {0}').format(fn))

    def test_versus_reference_results_split_pca_no_cov_normalization(self):
        self.run_command(tags=["--split_pca", "-i 100", "-m 1", "--no_cov_normalization"])
        reference_result_dir = os.path.join(test_dir_path,
                                            'test_data',
                                            'reference_result_split_pca_no_cov_normalization')
        for ref_f in listdir(reference_result_dir):
            fn = os.path.basename(ref_f)
            # Log will have time stamps and thus not the same
            if fn == 'log.txt':
                continue
            new_f = os.path.join(tmp_basename_dir,fn)
            ref_f = os.path.join(reference_result_dir,fn)
            with open(ref_f,'r') as ref_fh:
                ref_f = ref_fh.read()
            with open(new_f,'r') as new_fh:
                new_f = new_fh.read()
            assert_true(ref_f == new_f,
                        msg=('File not consistent with '
                             'reference file {0}').format(fn))

    def test_diagonal_covariance_matrix(self):
        self.run_command(tags=["--covariance_type diag","-i 100",'-m 1'])
        assert_equal(self.c,0,
                     msg = ("Command exited with nonzero status "
                            "when ran with diagonal cov matrix"))
        fn = 'bic.csv'
        ref_f = os.path.join(test_dir_path, 'test_data',
                             'reference_result', fn)
        new_f = os.path.join(tmp_basename_dir,fn)
        with open(ref_f,'r') as ref_fh:
            ref_f = ref_fh.read()
        with open(new_f,'r') as new_fh:
            new_f = new_fh.read()
        assert_true(new_f != ref_f,
                    msg=('Diagonal cov matrix clustering consistent with '
                         'reference clustering.'))

