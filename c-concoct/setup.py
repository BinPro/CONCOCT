from distutils.core import setup, Extension
 
module1 = Extension('vbgmm', 
	libraries =['gsl',  'gslcblas'],
	sources = ['vgmmmodule.c'])
 
setup (name = 'vbgmm',
        version = '1.0',
        description = 'write description here',
        ext_modules = [module1])
