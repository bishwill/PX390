#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>


// First we define all the functions and structures that we need from band_utlity.c

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


// Next we have functions to read the inputs and coefficients into arrays.

// Reads the parameters into the program from input.txt
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

// Reads the coefficients into the program from coefficients.txt
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

// This function will return the result of i mod N. This is needed for indexing when assigning values to the matrix.
// For example, for an array of length 10: mod(1, 10) = 1, mod(10, 10) = 0, mod(11, 10) = 1, mod(-1, 10) = 9
long mod(long i, long N) {
  return (i % N + N) % N;
}

// This function applies the matrix folding operation to an index i from an array of length N.
// This folding operation is described in Week 6 Lectures.
long fold_mapping(long i, double N) {
  if (i < N/2) {
    return 2 * i;
  }
  else {
    return 2 * (N - i) - 1;
  }
}

// Solution starts here
int main(void) {

  // Reading in input variables
  double L;
  long N;
  double D;
  double v;
  double kPlus;
  double kMinus;
  read_input(&L, &N, &D, &v, &kPlus, &kMinus);
  
  // Calculating the values required for the banded matrix
  long K = 2*N;
  double dx = L / N;
  double dx2 = dx*dx;
  double alpha = (D / dx2) + (v / dx);
  double beta = (-2*D / dx2) - (v / dx) - kPlus;
  double gamma = D / dx2;
  double delta = (-2*D / dx2) - (v / dx) - kMinus;


  // Reading coefficients into arrays. Note that S has been re-formulated to be length 2N.
  // Since it contains the values of S and zeros corresponding to the 2 equations.
  double *S = (double *) malloc(sizeof(double) * 2 * N);
  double *sigma = (double *) malloc(sizeof(double) * N);
  read_coefficients(S, sigma, N);
  
  // Applying the matrix fold mapping to S since the matrix is not narrow banded.
  double *mapped_S = (double *) malloc(sizeof(double) * 2 * N);
  for (long k = 0; k < N; k++) {
    mapped_S[fold_mapping(k, K)] = S[k];
  }
  

  // Initialising the band matrix
  // With this implementation we get 4 upper bands and 4 lower bands
  band_mat bmat;
  long nbands_low = 4;
  long nbands_up  = 4;
  init_band_mat(&bmat, nbands_low, nbands_up, K);

  
  // Setting the banded matrix using setv() from band_utility.c
  // The first step is to take the index i modulo K which will handle indices that are out of bounds.
  // Then I apply the fold_mapping function to turn the matrix into narrow banded form.
  // Then use this index in the setv() function to create the banded matrix.
  for (long i = 0; i < N; i++) {
    // Setting alpha
    setv(&bmat, fold_mapping(2*i, K), fold_mapping(mod(2*i-2, K), K), alpha);
    setv(&bmat, fold_mapping(2*i+1, K), fold_mapping(mod(2*i-1, K), K), alpha);
    
    // Setting beta
    setv(&bmat, fold_mapping(2*i, K), fold_mapping(2*i, K), beta);
    

    // Setting kPlus
    setv(&bmat, fold_mapping(2*i+1, K), fold_mapping(2*i, K), kPlus);

    // Setting kMinus
    setv(&bmat, fold_mapping(2*i, K), fold_mapping(2*i+1, K), kMinus);
    
    // Setting gamma
    setv(&bmat, fold_mapping(2*i, K), fold_mapping(mod(2*i+2, K), K), gamma);
    setv(&bmat, fold_mapping(2*i+1, K), fold_mapping(mod(2*i+3, K), K), gamma);

    // Setting delta
    setv(&bmat, fold_mapping(2*i+1, K), fold_mapping(2*i+1, K), delta - sigma[i]);
  }

  // Creating an array with the output values and solving the system.
  double *mapped_x = (double *) malloc(sizeof(double) * 2 * N);
  solve_Ax_eq_b(&bmat, mapped_x, mapped_S);


  // Applying the inverse mapping to get the output values in the correct order.
  double *x = (double *) malloc(sizeof(double) * 2 * N);
  for (long k = 0; k < K; k++) {
    x[k] = mapped_x[fold_mapping(k, K)];
  }

  // Writing x, A(x) and B(x) to output.txt
  FILE *output_file;
  output_file = fopen("output.txt","w");
  if(output_file == NULL) {
    printf("Error writing to file\n");
    exit(1);
  }
  for (long k = 0; k < N; k++) {
    fprintf(output_file, "%lf %lf %lf\n", k*dx, x[2*k], x[2*k + 1]);
  }
  fclose(output_file);

  
  // Freeing up used memory
  free(S);
  free(mapped_S);
  free(x);
  free(mapped_x);
  free(sigma);
  finalise_band_mat(&bmat);
  return 0;
}