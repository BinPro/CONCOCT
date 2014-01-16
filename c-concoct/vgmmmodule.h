#ifndef NMGS_H
#define NMGS_H

typedef struct s_Params
{
  /*seed*/
  unsigned long int lSeed;
  /*csv input file*/
  char *szInputFile;
  /*pca transformation*/
  char *szPInputFile;
  /*output file stub*/
  char *szOutFileStub;
  /*initial cluster size*/
  int nKStart;
  /*min contig length*/
  int nLMin;
} t_Params;


typedef struct s_Data
{
  int nN;

  int nT;

  int nD;

  gsl_matrix *ptTMatrix;

  double **aadX;

  char **aszDimNames;

  char **aszSampleNames;
} t_Data;

typedef struct s_VBParams
{
  /*scale for mean prior*/
  double dBeta0;
  
  /*Wishart degrees of freedom*/
  double dNu0;

  /*Inverse! of the Wishart scale precision-matrix*/
  gsl_matrix *ptInvW0;

  /*Log Wishart normalisation*/
  double dLogWishartB;

} t_VBParams;


typedef struct s_Cluster
{
  /*parameters for variational Bayes*/
  t_VBParams *ptVBParams;
  /*start seed*/
  unsigned long lSeed;
  /*thread index*/
  int nThread;
  /*pointer to data*/
  t_Data *ptData;
  /*number of data points*/
  int nN;
  /*number of clusters*/
  int nK;
  /*number of dimensions*/
  int nD;
  /*variational lower bound*/
  double dVBL;
  /*Means*/
  double **aadMu;
  /*Scaled weight Bishop 10.60*/
  double *adBeta;  
  /*Scaled means Bishop 10.61*/
  double **aadM;
  /*sample covariance matrix for each cluster storing this helps with lower bound calcn*/
  gsl_matrix **aptCovar;
  /*Inverse regularised variances Bishop 10.62*/
  gsl_matrix **aptSigma;
  /*Bishop 10.63*/
  double *adNu;
  /*Responsibilities*/
  double **aadZ;
  /*log-Matrix determinants*/
  double *adLDet;
  /*mixture weights*/
  double *adPi;
  /*assigned cluster for each data point*/
  int *anMaxZ;
  /*frequencies for each cluster*/
  int *anW;
} t_Cluster;


#define DELIM ",\n"
#define MAX_LINE_LENGTH   10000
#define MAX_WORD_LENGTH   128

#define TRUE  1
#define FALSE 0

#define NOT_SET -1
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

/*Default parameters*/
#define DEF_KSTART       400
#define DEF_LMIN         1000
#define DEF_BETA0        1.0e-3

#define MAX_ITER         500
#define MIN_CHANGE_VBL   1.0e-6
#define MIN_PI           0.1 /*Unormalised*/
#define MIN_COVAR        0.001

#define N_RTHREADS       10

#define K_PRIME          100003
#define R_PRIME          1009

#define DEF_SEED         1
#define DEF_LMIN         1000

#define OUT_FILE_STUB    "-out"
#define INPUT_FILE       "PCA_transformed_data_gt"
#define PINPUT_FILE      "PCA_components_data_gt"

int driver(const char* szFileStub);

void setParams(t_Params *ptParams,char *szFileStub);

void destroyParams(t_Params *ptParams);

void setVBParams(t_VBParams *ptVBParams, t_Data *ptData);

void readInputData(const char *szFile, t_Data *ptData);

void readPInputData(const char *szFile, t_Data *ptData);

void destroyData(t_Data *ptData);

void allocateCluster(t_Cluster *ptCluster, int nN, int nK, int nD, t_Data *ptData, long lSeed);

void performMStep(t_Cluster *ptCluster, t_Data *ptData);

void updateMeans(t_Cluster *ptCluster, t_Data *ptData);

void gmmTrainVB(t_Cluster *ptCluster, t_Data *ptData);

void initRandom(gsl_rng *ptGSLRNG, t_Cluster *ptCluster, t_Data *ptData);

void initKMeans(gsl_rng *ptGSLRNG, t_Cluster *ptCluster, t_Data *ptData);

double calcDist(double* adX, double *adMu, int nD);

void calcZ(t_Cluster* ptCluster, t_Data* ptData);

void writeClusters(char *szOutFile, t_Cluster *ptCluster, t_Data *ptData);

void destroyCluster(t_Cluster* ptCluster);

void* fitEM(void *pvCluster);

void* runRThreads(void *pvpDCluster);

void writeMeans(char *szOutFile, t_Cluster *ptCluster);

void writeTMeans(char *szOutFile, t_Cluster *ptCluster,t_Data *ptData);

void writeSquareMatrix(char*szFile, gsl_matrix *ptMatrix, int nD);

void calcSampleVar(t_Data *ptData,double *adVar, double *adMu);

double dLogWishartB(gsl_matrix *ptInvW, int nD, double dNu, int bInv);

double calcVBL(t_Cluster* ptCluster);

void compressCluster(t_Cluster* ptCluster);

void calcCovarMatrices(t_Cluster *ptCluster, t_Data *ptData);

#endif
