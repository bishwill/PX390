/*********************************************************
* This program solves coupled equations which describs
* the diffusion-damped system with constant coefficients.
* A time splitting technique by Marchuk–Strang is used. 
* For two subproblems, A1 and A2 we solve A1 with time length
* tau=(1/2)*dt then solve the other one, A2, for time length
* dt and solve with A1 again for time length tau=(1/2)*dt.
* 
* The read_input function should not be modified, but the
* use of the function in main() may or may not be correct.
*
* To compile: gcc -Wall -Werror -std=c99 -o assign3 assign3.c -lm
*
* List of identified errors:
*--------+--------------------------------------------
*  Line  |     Brief description of a fix
* Number |

    32         Changed the output of the main function
               from double to int because the main
               function must always return an output
               of type int

    51         Don't know what to put here

    98         Function sin() was called but it hasn't
               been defined due to the fact that the
               math library has not been included. To
               fix, I have added the line
               #include <math.h> to the program so sin()
               is defined.

    109        Function cos() was called but it hasn't
               been defined due to the fact that the
               math library has not been included. To
               fix, I have added the line
               #include <math.h> to the program so cos()
               is defined.

    96-97      Within the if statement, the expressions
               used = instead of ==. This is an
               assignment rather than a logical
               expression so it should not be used in an
               if statement. Changing each = to a == will
               fix the issue.

    110        ctime was declared but it was not given
               an initial value. Therefore when it was
               incremented further in the program, it
               did not have a value to increment. To fix
               it, I set the value to be 0.0 when
               the variable is declared.

*--------+-------------------------------------------
* Example (not a real error):
*  21 ...... Removed space between(void) and {
********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI  3.141592

void read_input(double *D, double *L, int *nx, double *t_F);

int main(void) {
  /******************************/
  /* Declarations of parameters */
  /******************************/
  /* Number of grid points */                  
  int nx;
  /* Length of domain */
  double L;
  /* Equation coefficients */
  double D;
  /* Length of time to run simulation. */
  double t_F;

  /* Read in from file; */
  read_input(&D, &L, &nx, &t_F);

  /* Grid spacing */
  double dx = L/nx;
  double invdx2 = 1.0/(dx*dx);      
  /* Time step */
  double dt = 0.5;

  /************************************************/
  /* Solution Storage at Current / Next time step */
  /************************************************/
  double *uc, *un, *vc, *vn;
  /* Time splitting solutions */
  double *uts1, *uts2, *vts1, *vts2;
  /* Derivative used in finite difference */
  double deriv;

  /* Allocate memory according to size of nx */
  uc = malloc(nx);
  un = malloc(nx);
  vc = malloc(nx);
  vn = malloc(nx);
  uts1 = malloc(nx);
  uts2 = malloc(nx);
  vts1 = malloc(nx);
  vts2 = malloc(nx);

  /* Check the allocation pointers */
  if (uc==NULL||un==NULL||vc==NULL||vn==NULL||uts1==NULL||
  uts2==NULL||vts1==NULL||vts2==NULL) {
    printf("Memory allocation failed\n");
    return 1;
  }
  
  int k;
  double x;
  /* Current time */
  double ctime = 0.0;

  /* Initialise arrays */
  for(k = 0; k <= nx; k++) {
    x = k*dx;
    uc[k]  = 1.0 + sin(2.0*PI*x/L);
    vc[k]  = 0.0;
    /* Set other arrays to 0 */
    uts1[k] = 0; uts2[k] = 0;
    vts1[k] = 0; vts2[k] = 0;
  }

  /* Loop over timesteps */ 
  while (ctime < t_F){
    
    /* Rotation factors for time-splitting scheme. */
    double cfac = cos(2*dt);
    double sfac = sin(2*dt);
    
    /* First substep for diffusion equation, A_1 */ 
    for (k = 0; k < nx; k++) {
      x = k*dx;
      /* Diffusion at half time step. */
      deriv = (uc[k-1] + uc[k+1] - 2*uc[k])*invdx2;
      uts1[k] = uc[k] + D * deriv * 0.5*dt;
      deriv = (vc[k-1] + vc[k+1] - 2*vc[k])*invdx2;
      vts1[k] = vc[k] + D * deriv * 0.5*dt;
    }

    /* Second substep for decay/growth terms, A_2 */
    for (k = 0; k < nx; k++) {
      x = k*dx;
      /* Apply rotation matrix to u and v, */
      uts2[k] = cfac*uts1[k] + sfac*vts1[k];
      vts2[k] = sfac*uts1[k] + cfac*vts1[k];
    }

    /* Thirs substep for diffusion terms, A_1 */
    for (k = 0; k < nx; k++) {
      x = k*dx;
      deriv = (uts2[k-1] + uts2[k+1] - 2*uts2[k])*invdx2;
      un[k] = uts2[k] + D * deriv * 0.5*dt;
      deriv = (vts2[k-1] + vts2[k+1] - 2*vts2[k])*invdx2;
      vn[k] = vts2[k] + D * deriv * 0.5*dt;
    }
	   
    /* Copy next values at timestep to u, v arrays. */
    double *tmp;
    tmp = vts1;
    vn = vts2; 
    vc = tmp;
    tmp = un;
    un = uc; 
    uc = tmp;

    /* Increment time. */   
    ctime += dt;
    for (k = 0; k < nx; k++ ) {
      x = k*dx;
      printf("%g %g %g %g\n",ctime,x,uc[k],vc[k]);
    }
  }
  /* Free allocated memory */
  free(&uc); free(&un);
  free(&vc); free(&vn);
  free(&uts1); free(&uts2);
  free(&vts1); free(&vts2);

  return 0;
}

// The lines below don't contain any bugs! Don't modify them
void read_input(double *D, double *L, int *nx, double *t_F) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(4!=fscanf(infile,"%lf %lf %d %lf",D,L,nx,t_F)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}