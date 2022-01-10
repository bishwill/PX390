#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>

/* Define structure that holds band matrix information */
struct band_mat{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
/* Define a new type band_mat */
typedef struct band_mat band_mat;

int solve_Ax_eq_b(band_mat *bmat, double *x, double *b);
int printmat(band_mat *bmat);
long mod(long i, long N);
void read_input(double *L, long *N, double *D, double *v, double *kPlus, double *kMinus);
void read_coefficients(double *S, double *sigma, long N);
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns);
void finalise_band_mat(band_mat *bmat);
double *getp(band_mat *bmat, long row, long column);
double getv(band_mat *bmat, long row, long column);
double setv(band_mat *bmat, long row, long column, double val);
long fold_mapping(long i, double N);
long fold_mapping_inv(long i, double N);

int main(void) {

  // Input values
  double L;
  long N;
  double D;
  double v;
  double kPlus;
  double kMinus;
  read_input(&L, &N, &D, &v, &kPlus, &kMinus);
  
  double dx = L / N;
  double dx2 = dx*dx;
  double alpha = (D / dx2) + (v / dx);
  double beta = (-2*D / dx2) - (v / dx) - kPlus;
  double gamma = D / dx2;
  double delta = (-2*D / dx2) - (v / dx) - kMinus;


  // Coefficients
  double *S = (double *) malloc(sizeof(double) * 2 * N);
  double *sigma = (double *) malloc(sizeof(double) * N);
  read_coefficients(S, sigma, N);
  
  double *mapped_S = (double *) malloc(sizeof(double) * 2 * N);
  for (long k = 0; k < N; k++) {
    mapped_S[fold_mapping(k, 2*N)] = S[k];
  }
  

  // // Printing out values
  // printf("Input Values:\n");
  // printf("L = %lf\n", L);
  // printf("N = %ld\n", N);
  // printf("D = %lf\n", D);
  // printf("v = %lf\n", v);
  // printf("kPlus = %lf\n", kPlus);
  // printf("kMinus = %lf\n", kMinus);
  // printf("dx = %lf\n", dx);

  // printf("alpha = %lf\n", alpha);
  // printf("beta = %lf\n", beta);
  // printf("gamma = %lf\n", gamma);
  // printf("delta = %lf\n", delta);

  // printf("\n");

  // printf("Coefficients:\n");
  // for (int i = 0; i < 2*N; i++) {
  //   printf("S_%d = %lf\n", i, S[i]);
  // }
  // for (int i = 0; i < N; i++) {
  //   printf("sigma_%d = %lf\n", i, sigma[i]);
  // }


  

  band_mat bmat;
  // With this implementation we get 4 upper bands and 4 lower bands
  long nbands_low = 4;
  long nbands_up  = 4;
  
  init_band_mat(&bmat, nbands_low, nbands_up, 2*N);

  

  for (long i = 0; i < N; i++) {
    
    // printf("%ld %ld %lf\n", 2*i, mod(2*i-2, 2*N), alpha);
    // printf("%ld %ld %lf\n", 2*i+1, mod(2*i-1, 2*N), alpha);

    // printf("%ld %ld %lf\n", 2*i, mod(2*i, 2*N), beta);

    // printf("%ld %ld %lf\n", 2*i+1, mod(2*i, 2*N), kPlus);

    // printf("%ld %ld %lf\n", 2*i, mod(2*i+1, 2*N), kMinus);

    // printf("%ld %ld %lf\n", 2*i, mod(2*i+2, 2*N), gamma);
    // printf("%ld %ld %lf\n", 2*i+1, mod(2*i+3, 2*N), gamma);

    // printf("%ld %ld %lf\n", 2*i+1, mod(2*i+1, 2*N), delta - sigma[i]);




    // printf("%ld %ld %lf\n", fold_mapping(2*i, 2*N), fold_mapping(mod(2*i-2, 2*N), 2*N), alpha);
    // printf("%ld %ld %lf\n", fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i-1, 2*N), 2*N), alpha);

    // printf("%ld %ld %lf\n", fold_mapping(2*i, 2*N), fold_mapping(mod(2*i, 2*N), 2*N), beta);

    // printf("%ld %ld %lf\n", fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i, 2*N), 2*N), kPlus);

    // printf("%ld %ld %lf\n", fold_mapping(2*i, 2*N), fold_mapping(mod(2*i+1, 2*N), 2*N), kMinus);

    // printf("%ld %ld %lf\n", fold_mapping(2*i, 2*N), fold_mapping(mod(2*i+2, 2*N), 2*N), gamma);
    // printf("%ld %ld %lf\n", fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i+3, 2*N), 2*N), gamma);

    // printf("%ld %ld %lf\n", fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i+1, 2*N), 2*N), delta - sigma[i]);


    // Setting alpha
    setv(&bmat, fold_mapping(2*i, 2*N), fold_mapping(mod(2*i-2, 2*N), 2*N), alpha);
    setv(&bmat, fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i-1, 2*N), 2*N), alpha);
    
    // Setting beta
    setv(&bmat, fold_mapping(2*i, 2*N), fold_mapping(mod(2*i, 2*N), 2*N), beta);
    

    // Setting kPlus
    setv(&bmat, fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i, 2*N), 2*N), kPlus);

    // Setting kMinus
    setv(&bmat, fold_mapping(2*i, 2*N), fold_mapping(mod(2*i+1, 2*N), 2*N), kMinus);
    
    // Setting gamma
    setv(&bmat, fold_mapping(2*i, 2*N), fold_mapping(mod(2*i+2, 2*N), 2*N), gamma);
    setv(&bmat, fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i+3, 2*N), 2*N), gamma);

    // Setting delta
    setv(&bmat, fold_mapping(2*i+1, 2*N), fold_mapping(mod(2*i+1, 2*N), 2*N), delta - sigma[i]);
  }


  // printmat(&bmat);

  double *mapped_x = (double *) malloc(sizeof(double) * 2 * N);
  solve_Ax_eq_b(&bmat, mapped_x, mapped_S);

  // for (long k = 0; k < N; k++) {
  //   printf("(%lf, %lf)\n", k*dx, mapped_x[2*k]);
  // }

  double *x = (double *) malloc(sizeof(double) * 2 * N);
  for (long k = 0; k < 2*N; k++) {
    x[fold_mapping_inv(k, 2*N)] = mapped_x[k];
  }
  

  for (long k = 0; k < N; k++) {
    printf("%lf %lf %lf\n", k*dx, x[2*k], x[2*k + 1]);
  }


  free(S);
  free(mapped_S);
  free(x);
  free(sigma);
  finalise_band_mat(&bmat);
  return 0;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) { 
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
      
    }
  }
  return 0;
}

long mod(long i, long N) {
  return (i % N + N) % N;
}

void read_input(double *L, long *N, double *D, double *v, double *kPlus, double *kMinus) {
  FILE *infile;
  if(!(infile=fopen("input.txt","r"))) {
    printf("Error opening file\n");
    exit(1);
  }
  if(6 != fscanf(infile, "%lf %ld %lf %lf %lf %lf", L, N, D, v, kPlus, kMinus)) {
    printf("Error reading parameters from file\n");
    exit(1);
  }
  fclose(infile);
}

void read_coefficients(double *S, double *sigma, long N) {
  FILE *coeff_file;
  if(!(coeff_file=fopen("coefficients.txt","r"))) {
    printf("Error opening file\n");
    exit(1);
  }
  for(long j = 0; j < N; j++) {
    if (2 != fscanf(coeff_file, "%lf %lf", &S[2*j], &sigma[j])) {
      printf("Error reading coefficients from file\n");
      exit(1);
    }
    S[2*j] = -S[2*j];
    S[2*j+1] = 0;
  }
  fclose(coeff_file);
}

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */ 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }  
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
}

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
}

long fold_mapping(long i, double N) {
  if (i < N/2) {
    return 2 * i;
  }
  else {
    return 2 * (N - i) - 1;
  }
}

long fold_mapping_inv(long i, double N) {
  if (i % 2 == 0) {
    return i/2;
  }
  else {
    return (-i - 1 + 2*N) / 2;
  }
}

