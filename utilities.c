#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_vector(int nrow, double *a)
{
  int i;
  for(i=0; i<nrow; i++)
    printf(" %.6e ", a[i]);
  printf("\n");
}

void create_matrix(int nxglob, double **a, double gamma)
{
	int i,j;

	for(i=0; i<nxglob; i++)
		for(j=0; j<nxglob; j++)
			a[i][j] = 0;
		
	for(i=1; i<nxglob-1; i++)
	{
		a[i][i-1] = -gamma;	
		a[i][i] = 1+2*gamma;
		a[i][i+1] = -gamma;
	}
	a[0][0] = 1;
	a[nxglob-1][nxglob-1] = 1;
}	
	

void grid(int nx, int nxglob, int istglob, int ienglob, double xstglob, double xenglob, double *x, double dx)
{
  int i, iglob;

 // *dx = (xenglob - xstglob)/(double)(nxglob-1);

  for(i=0; i<nx; i++)
  {
    iglob = istglob + i;
    x[i] = xstglob + (double)iglob * (dx);
  }
}

void enforce_bcs(int nx, int nxglob, int istglob, int ienglob, double *x, double *T)
{
  if(istglob==0)
    T[0] = -1.0;

  if(ienglob==nxglob)
  T[nx-1] = 1.0;
}

void get_exact_solution(int nx, double time, double *x, double *Texact)
{
  int i;

  for(i=0; i<nx; i++)
    Texact[i] = erf((x[i]-0.5)/2.0/sqrt(time));
}

void set_initial_condition(int nx, int nxglob, int istglob, int ienglob, double time, double *x, double *T)
{
  int i;

  get_exact_solution(nx, time, x, T);

  enforce_bcs(nx,nxglob,istglob,ienglob,x,T);
}

void vec_norm(int nx, int nxloc, double *bp, double *locnormb)
{
  int i;
 // double locnormb = 0.0;
   *locnormb =  0.0; 
  for(i=0; i<nxloc; i++)
    *locnormb += bp[i]*bp[i];

// MPI_Allreduce(&locnormb, &normb,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

//  *normb = sqrt(*normb/(double)nx);
}

void dot_prod(int nxloc, int istglob, int ienglob, double *v1, double *v2, double *locresult)
{
  int i;
  *locresult = 0.0;
  for(i=0; i<nxloc; i++)
    *locresult += v1[i]*v2[i];
 //  *result = locresult;
 //  MPI_Allreduce(&locresult, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void matrix_multiply(int nx, int nxloc, int istglob, int ienglob, double **ap, double *xp, double *bp, double *x)
{
  int i, j,iloc;
 // double *x;
 // x = (double *)malloc(nx*sizeof(double));
//  MPI_Allgather(xp,nxloc,MPI_DOUBLE,x,nxloc,MPI_DOUBLE,MPI_COMM_WORLD);
  for(i=istglob; i<ienglob; i++)
  {
    iloc = i-istglob;
    bp[iloc] = 0.0;
    for(j=0; j<nx; j++)
      bp[iloc] += ap[iloc][j]*x[j];
  }
}

void output_soln(int rank, int nx, int it, double tcurr, double *x, double *T, double *Tex)
{
  int i;
  FILE* fp;
  char fname[100];

  sprintf(fname, "T_x_%04d_%04d.dat", it, rank);

  fp = fopen(fname, "w");
  for(i=0; i<nx; i++)
    fprintf(fp, "%.15e %.15e %.15e\n", x[i], T[i], Tex[i]);
  fclose(fp);
}