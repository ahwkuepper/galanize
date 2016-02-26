#define Pi 3.14159265
#define PI 3.14159265
#define G  0.0043           //gravitational constant in [km^2/s^2/Msun*pc]
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define NMAX 10000010  //maximum allowed particle number
#define ROTATE(am,i,j,k,l) g=am[i][j];h=am[k][l];am[i][j]=g-s*(h+g*tau);am[k][l]=h+s*(g-h*tau);
#define NR_END 1
#define FREE_ARG char*

void radial_profile(double **star, int N, double M, double Rscale, double *Rh2D, double *Rh3D);
void triax(double **particle, int NLAGR, double *ha, double *hb, double *hc, double *LAGR, int nbody);
void beta(double **particle, double Rscale, int nbody);
void bootstrap(double **aa, int NR, double *er1, double *er2);
void shellsort_reverse_1d(double *array, int N);
void shellsort_1d(double *array, int N);
void shellsort(double **array, int N, int k);
void shellsort_reverse(double **array, int N, int k);

void jacobi(float **am, int n, float d[], float **v, int *nrot);
void eigsrt(float d[], float **v, int n);
void nrerror(char error_text[]);
float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);

