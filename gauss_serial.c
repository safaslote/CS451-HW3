/* Gaussian elimination without pivoting.
 * Compile with "gcc gauss_serial.c" 
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#include <stdio.h>

/* Program Parameters */
#define MAXN 2000  /* Max value of N */
int N;  /* Matrix size */
int numThreads; /* Number of threads */
int npes; //number of processors
int myrank; //rank
int start;// start for wtime function
int end;//end for wtime function

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */

/* junk */
// #define randm() 4|2[uid]&3

/* Prototype */
void gauss();  /* The function you will provide.
		* It is this routine that is timed.
		* It is called only on the parent.
		*/

/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
	int seed = 0;  /* Random seed */
	// char uid[32]; /*User name */

	/* Read command-line arguments */
	srand(time_seed());  /* Randomize */

	if (argc == 4) {
		seed = atoi(argv[2]);
		srand(seed);
		numThreads=atoi(argv[3]);
		printf("Random seed = %i\n", seed);
	} 
	if (argc >= 2) {
		N = atoi(argv[1]);
		if (N < 1 || N > MAXN) {
			printf("N = %i is out of range.\n", N);
			exit(0);
		}
	}
	else {
		printf("Usage: %s <matrix_dimension> [random seed]\n",
				argv[0]);    
		exit(0);
	}

	/* Print parameters */
	printf("\nMatrix dimension N = %i.\n", N);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}

/* Print input matrices */
void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
	printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

void print_X() {
  int row;

  if (N < 100) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}

int main(int argc, char **argv) {
  int i;
  /* Timing variables */
  struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
  struct timezone tzdummy;
  // clock_t etstart2, etstop2;  /* Elapsed times using times() */
  unsigned long long usecstart, usecstop;
  struct tms cputstart, cputstop;  /* CPU times for my processes */

  //initialize MPI
  MPI_Init(&argc, &argv);
  //Obtain partition size
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  //Obtain number of processors
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  /* Process program parameters */
  parameters(argc, argv);
  //data needs to be held and generated in process 0, so generate inputs in root 0
  if(npes == 0){
	initialize_inputs();
  }
  //where we do gaussian elimination calculations
  gauss();
  
 
  MPI_Finalize(); //end of MPI computation
  /* Start Clock */
  printf("\nStarting clock.\n");
  gettimeofday(&etstart, &tzdummy);
  times(&cputstart);


  /* Stop Clock */
  gettimeofday(&etstop, &tzdummy);
  times(&cputstop);
  printf("Stopped clock.\n");
  usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
  usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

  /* Display output */
  print_X();

  /* Display timing results */
  printf("\nElapsed time = %g ms.\n",
	 (float)(usecstop - usecstart)/(float)1000);

  printf("(CPU times are accurate to the nearest %g ms)\n",
	 1.0/(float)CLOCKS_PER_SEC * 1000.0);
  printf("My total CPU time for parent = %g ms.\n",
	 (float)( (cputstop.tms_utime + cputstop.tms_stime) -
		  (cputstart.tms_utime + cputstart.tms_stime) ) /
	 (float)CLOCKS_PER_SEC * 1000);
  printf("My system CPU time for parent = %g ms.\n",
	 (float)(cputstop.tms_stime - cputstart.tms_stime) /
	 (float)CLOCKS_PER_SEC * 1000);
  printf("My total CPU time for child processes = %g ms.\n",
	 (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
		  (cputstart.tms_cutime + cputstart.tms_cstime) ) /
	 (float)CLOCKS_PER_SEC * 1000);
      /* Contrary to the man pages, this appears not to include the parent */
  printf("--------------------------------------------\n");
  
  exit(0);
}

/* ------------------ Above Was Provided --------------------- */


/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */
void gauss() {
  int norm, row, col;  /* Normalization row, and zeroing
			* element row and col */
  float multiplier;
MPI_Status status; //need status in order to use MPI receive function
MPI_Request request;

//if root = 0, start the time since that is where we start the computation
if(myrank == 0){
	start = MPI_Wtime();
	printf("MPI Start Time %d\n", start);
}

  /* Gaussian elimination */
  for (norm = 0; norm < N - 1; norm++) {
 //we need to send the normalized data for each iteration 
 //norm row should be applied to all processors for each row of the matrix
 //do it for both the matrix and the vector
  MPI_Bcast(&A[norm], N, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&B[norm], 1, MPI_INT, 0, MPI_COMM_WORLD); 
   //use MPI scatter to assign each row to a processor
    for (row = norm + 1; row < N; row++) {
       multiplier = A[row][norm] / A[norm][norm];
      for (col = norm; col < N; col++) {
	      A[row][col] -= A[norm][col] * multiplier;
      }
      B[row] -= B[norm] * multiplier;
     
 }
} 

  /* Back substitution */
  for (row = N - 1; row >= 0; row--) {
    X[row] = B[row];
    for (col = N-1; col > row; col--) {
      X[row] -= A[row][col] * X[col];
    }
    X[row] /= A[row][row];
  }
//if we are at root processor, end the Wtime since all data gets returned to the root
if(myrank == 0){
 end = MPI_Wtime();
 printf("End time %d\n",end);
}
}
