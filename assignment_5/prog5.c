// This program will solve the problem in the following way
// At each time step, it will implicitly solve the spatial derivative
// by solving a matrix equation of u_{t+1} in terms of u_{t}/
// and then print out the values of u at all grid points.
// It then increments the time step and repeats.

// Moodules required
#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>



// First we define the functions needed for this program from band_utility.c 

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

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  // printf("%ld %ld %lf\n", row, column, val);
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

// This function reads the coefficients in from coefficients.txt
void read_coefficients(double *u0, long N) {
  FILE *coeff_file;
  if(!(coeff_file=fopen("coefficients.txt","r"))) {
    printf("Error opening file\n");
    exit(1);
  }
  for(long k = 0; k < N; k++) {
    if (1 != fscanf(coeff_file, "%lf", &u0[k])) {
      printf("Error reading coefficients from file\n");
      exit(1);
    }
  }
  fclose(coeff_file);
}

// This function reads the input values from input.txt
void read_input(long *Nx, long *Ny, double *Lx, double *Ly, double *tf, double *tD, int *sx0, int *sx1, int *sy0, int *sy1) {
    FILE *infile;
    if(!(infile=fopen("input.txt","r"))) {
        printf("Error opening file\n");
        exit(1);
    }

    if (1 != fscanf(infile, "%ld", Nx)) {
        printf("Error reading Nx from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%ld", Ny)) {
        printf("Error reading Ny from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%lf", Lx)) {
        printf("Error reading Lx from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%lf", Ly)) {
        printf("Error reading Ly from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%lf", tf)) {
        printf("Error reading tf from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%lf", tD)) {
        printf("Error reading tD from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%d", sx0)) {
        printf("Error reading sx0 from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%d", sx1)) {
        printf("Error reading sx1 from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%d", sy0)) {
        printf("Error reading sy0 from file\n");
        exit(1);
    }
    if (1 != fscanf(infile, "%d", sy1)) {
        printf("Error reading sy1 from file\n");
        exit(1);
    }
  fclose(infile);
}

// This function gets the index value of the array that stores the grid points
// from a coordinate (i, j). It is used when printing out the values to output.txt
double get_array_index(long i, long j, long Nx, long Ny) {
  // Currently only does Dirichlet
  if (i >= Nx || i < 0) {
      printf("Row index out of range %ld\n", i);
      exit(1);
    }
    if (j >= Ny || j < 0) {
      printf("Column index out of range %ld\n", j);
      exit(1);
    }
  return j * Nx + i;
}

// This function does the inverse of the above function. It takes an array
// index value k and returns the x coordinate of the grid i.
long get_i(long k, long Nx) {
    return k % Nx;
}

// And this function will get the y coordinate of the grid given an array
// index value k
long get_j(long k, long Nx) {
    long i =  get_i(k, Nx);
    return (k-i) / Nx;
}

// Program starts here
int main(void) {

  // Variables to be read in from input.txt
  long Nx;
  long Ny;
  double Lx;
  double Ly;
  double tf;
  double tD;
  int sx0;
  int sx1;
  int sy0;
  int sy1;

  // Reading data from input.txt file into the variables
  read_input(&Nx, &Ny, &Lx, &Ly, &tf, &tD, &sx0, &sx1, &sy0, &sy1);
  
  // Allocating Memory

  // Grid point values of u_{t}
  double* u_current = (double *) malloc(sizeof(double) * Nx*Ny);
  // This will be the output array that solves the equation Ax = b.
  // It will contain the grid point values of u{t+1}
  double* x = (double *) malloc(sizeof(double) * Nx*Ny);
  // This is the b part of the system of equations Ax = b
  // It will be filled later on
  double* b = (double *) malloc(sizeof(double) * Nx*Ny);
  // Checking pointers have been allocated correctly
  if(u_current == NULL || x == NULL || b == NULL) {
    printf("Memory allocation failed\n");
    exit(1);
    }
  
  // Reading data from coefficients.txt into the array u_current
  read_coefficients(u_current, Nx*Ny);
  

  // Defining common values
  double dx = Lx / (Nx-1);
  double dy = Ly / (Ny-1);
  double invdx2 = 1 / (dx*dx);
  double invdy2 = 1 / (dy*dy);
  double invtD = 1 / tD;
  

  // Setting up matrix A
  // There will be Nx lower bands and Nx upper bands
  band_mat bmat;
  long nbands_low = Nx;
  long nbands_up  = Nx;
  init_band_mat(&bmat, nbands_low, nbands_up, Nx*Ny);

  // Setting values to A dependent on the boundary switches
  for (long k = 0; k < Nx*Ny; k++) {
    
    // Setting the main diagnoal (does not depend on switches)
    setv(&bmat, k, k, invtD + (2 * invdx2) + (2 * invdy2));

    // Setting the interior points (This will set the matrix for dirichlet conditions)
    if (get_i(k, Nx) > 0){
      setv(&bmat, k, k-1, -invdx2);
    }
    if (get_i(k, Nx) < Nx-1){
      setv(&bmat, k, k+1, -invdx2);
    }
    if (get_j(k, Nx) > 0){
      setv(&bmat, k, k-Nx, -invdy2);
    }
    if (get_j(k, Nx) < Ny-1){
      setv(&bmat, k, k+Nx, -invdy2);
    }
    

    // These 4 if statements will change certain matrix values depending
    // on the input switches (Van Neumann conditions)

    // Checking x = 0 boundary
    if (get_i(k, Nx) == 0 && sx0 == 1) {
      // u_{-1}{j} == u_{1}{j}
      setv(&bmat, k, k+1, -2*invdx2);
    }

    // Checking x = Lx boundary
    if (get_i(k, Nx) == Nx-1 && sx1 == 1) {
      // u_{Nx}{j} == u_{Nx-2}{j}
      setv(&bmat, k, k-1, -2*invdx2);
    }

    // Checking y = 0 boundary
    if (get_j(k, Nx) == 0 && sy0 == 1) {
      // u_{i}{-1} == u_{i}{1}
      setv(&bmat, k, k+Nx, -2*invdy2);
    }

    // Checking y = Ny boundary
    if (get_j(k, Nx) == Ny-1 && sy1 == 1) {
      // u_{i}{Nx} == u_{i}{Nx-2}
      setv(&bmat, k, k-Nx, -2*invdy2);
    }    
  }

  // Initialisiing time
  double ctime = 0.0;

  // Write initial values to output.txt
  FILE *output_file;
  output_file = fopen("output.txt","w");
  if(output_file == NULL) {
    printf("Error writing to file\n");
    exit(1);
  }
  for (long j = 0; j < Ny; j++) {
    for (long i = 0; i < Nx; i++) {
      long k = get_array_index(i, j, Nx, Ny);
      // printf("%ld\n", k);
      fprintf(output_file, "%lf %lf %lf %lf\n", ctime, i*dx, j*dy, u_current[k]);
    }
  }
  fclose(output_file);


  // Time evolution starts here
  // The while condition ensures ctime does not equal tf.
  // (For some values of tD, ctime would be equal to tf and it would enter the while
  //  loop when it shouldn't have so this condition prevents that from happening)
  while (ctime - tf < -0.0000001) {
  
    // Solve implicit matrix equation

    // Setting up b
    for (long k = 0; k < Nx*Ny; k++) {
      double u_k = u_current[k];
      double f_u_k = (u_k * u_k) / (1 + (u_k * u_k));
      b[k] = (u_k / tD) - f_u_k;
    }

    // Solving the system and outputting to x
    solve_Ax_eq_b(&bmat, x, b);

    // Increment time
    ctime += tD;

    // Writing grid point values to output.txt at time t+1
    FILE *output_file;
    output_file = fopen("output.txt","a");
    if(output_file == NULL) {
      printf("Error writing to file\n");
      exit(1);
    }
    for (long j = 0; j < Ny; j++) {
      for (long i = 0; i < Nx; i++) {
        long k = get_array_index(i, j, Nx, Ny);
        fprintf(output_file, "%lf %lf %lf %lf\n", ctime, i*dx, j*dy, x[k]);
        u_current[k] = x[k];
      }
    }
    fclose(output_file);
  }

  // Freeing memory and finalising band matrix and then exiting
  free(u_current);
  free(x);
  free(b);
  finalise_band_mat(&bmat);
  return 0;
}
