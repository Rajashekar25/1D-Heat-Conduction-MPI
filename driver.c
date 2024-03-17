#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void err_norm(int nx, int nxloc, double *x, double *xexact, int norm, double *normval)
{

  int i;
  double locnormval = 0.0;

  if(norm==1)
  {
    for(i=0; i<nxloc; i++)
      *normval += fabs(x[i]-xexact[i]);
    *normval /= (double) nx;
  }
  if(norm==2)
  {
    for(i=0; i<nx; i++)
      *normval += pow(x[i]-xexact[i], 2.0);
//  MPI_Allreduce(&locnormval,&normval,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    *normval = sqrt(*normval/(double) nx);
  }

}

int main(int argc, char** argv)
{

  int nx;
  double **a, **ap, *b, *bp, *workarrp; 
  double *x, *T, *Tp, *Texact, *Texactp;
  double tst, ten, xst, xen, dx, dt, tcurr, xlen, t_print, gamma, l2norm;
  int i, j, iloc, it, num_time_steps, it_print;
  char debugfname[100],time[100];
  FILE* fid, *fp, *fp2, *fp3;

  int rank, size;
  int nxglob, istglob, ienglob;
  double xstglob, xenglob;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(rank==0)
  {
    fp = fopen("input.in", "r");
    fscanf(fp, "%d\n", &nxglob);
    fscanf(fp, "%lf %lf\n", &xstglob, &xenglob);
    fscanf(fp, "%lf %lf\n", &tst, &ten);
    fscanf(fp, "%lf %lf\n", &dt, &t_print);
    fclose(fp);

    nx = nxglob/size;
    xlen = (xenglob-xstglob)/(double)size;
   dx = (xenglob - xstglob)/(double)(nxglob-1);

    num_time_steps = (int)((ten-tst)/dt) + 1;
    it_print = (int) (t_print/dt);
	gamma = dt/(dx*dx);
	
	  a = (double **) malloc(nxglob*sizeof(double *));
    for(i=0; i<nxglob; i++)
    a[i] = (double *) malloc(nxglob*sizeof(double));
  
   create_matrix(nxglob, a ,gamma);
   
   
   char vec[100],mat[100];  
   sprintf(mat, "input_matrices.in");
   fid = fopen(mat, "w");
//	fprintf(fid,"A matrix \n");
    for(i=0; i<nxglob; i++)
    { for(j=0; j<nxglob; j++)		
      fprintf(fid, "%lf ", a[i][j]);
      fprintf(fid,"\n");
	}
  fclose(fid);

  }


  int *sendarr_int;
  sendarr_int = malloc(4*sizeof(int));
  if(rank==0)
  {
    sendarr_int[0] = nxglob;         sendarr_int[1] = nx;
    sendarr_int[2] = num_time_steps; sendarr_int[3] = it_print;
  }
  MPI_Bcast(sendarr_int, 4, MPI_INT, 0, MPI_COMM_WORLD);
  if(rank!=0)
  {
            nxglob = sendarr_int[0];       nx = sendarr_int[1];
    num_time_steps = sendarr_int[2]; it_print = sendarr_int[3];
  }
  free(sendarr_int);

  double *sendarr_dbl;
  sendarr_dbl = malloc(9*sizeof(double));
  if(rank==0)
  {
    sendarr_dbl[0] = tst;  sendarr_dbl[1] = ten;     sendarr_dbl[2] = dt;      sendarr_dbl[3] = t_print;
    sendarr_dbl[4] = xlen; sendarr_dbl[5] = xstglob; sendarr_dbl[6] = xenglob; sendarr_dbl[7] = dx;
	sendarr_dbl[8] = gamma;
  }
  MPI_Bcast(sendarr_dbl, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(rank!=0)
  {
     tst = sendarr_dbl[0];     ten = sendarr_dbl[1];      dt = sendarr_dbl[2];  t_print = sendarr_dbl[4];
    xlen = sendarr_dbl[4]; xstglob = sendarr_dbl[5]; xenglob = sendarr_dbl[6]; dx = sendarr_dbl[7] ;
    gamma = sendarr_dbl[8];
  }
  free(sendarr_dbl);

 //  printf("dx-%d--%lf ",rank,dx);

  istglob = rank * (nxglob/size);
  ienglob = (rank+1) * (nxglob/size);

  xst = xstglob + rank*xlen;
  xen = xst + xlen;

  x = (double *)malloc(nx*sizeof(double));
  Tp = (double *)malloc(nx*sizeof(double));
//  rhs = (double *)malloc(nx*sizeof(double));
  Texactp = (double *)malloc(nx*sizeof(double));

  ap = (double **) malloc(nx*sizeof(double *));    //processor specific A matrix
  for(i=0; i<nx; i++)
    ap[i] = (double *) malloc(nxglob*sizeof(double));

    bp = (double *) malloc(nx*sizeof(double));
    workarrp = (double *) malloc(3*nx*sizeof(double));
    Texact = (double *)malloc(nxglob*sizeof(double));
    T = (double *)malloc(nxglob*sizeof(double));

  double dummyvar;
  
  fp3 = fopen("input_matrices.in","r");
     for(i=0; i<nxglob; i++)
     {
      if(i>=istglob && i<ienglob)
         {iloc = i-istglob;
          for(j=0;j<nxglob;j++)
	      fscanf(fp3,"%lf ",&ap[iloc][j]);
         }
      else
      {   for(j=0;j<nxglob;j++)
		fscanf(fp3,"%lf ",&dummyvar);
      }
    }
   fclose(fp3);

  
  grid(nx,nxglob,istglob,ienglob,xstglob,xenglob,x,dx);

  set_initial_condition(nx,nxglob,istglob,ienglob,tst,x,Tp);

  get_exact_solution(nx,tst,x,Texactp);

  output_soln(rank,nx,0,tst,x,Tp,Texactp);

 
  

  sprintf(debugfname, "debug_%04d.dat", rank);
  fid = fopen(debugfname, "w");
  fprintf(fid, "\n\n\n--Debug-1- %d %d %d %d %d\n", rank, nx, nxglob, istglob, ienglob);
  fprintf(fid, "\n-numtimesteps -%d it_print %d\n",num_time_steps,it_print);
  fprintf(fid, "--Debug-2- %d %lf %lf %lf %lf %lf\n", rank, xst, xen, xstglob, xenglob, xlen);
  fprintf(fid, "\n--Writing grid points--\n");
  for(i=0; i<nx; i++)
  {
    fprintf(fid, "%d %d %d %lf\n", rank, i, i+istglob, x[i]);
  }
  fprintf(fid, "--Done writing grid points--\n");
  fclose(fid);

  // initialize guess solution
 //  for(i=0; i<nx; i++)
 //    Tp[i] = 0.1;

  for(it=0; it<num_time_steps; it++)
  {
    tcurr = tst + (double)(it) * dt;

   for(i=0; i<nx; i++)  
    bp[i] = Tp[i];

    conjug_grad(nxglob, nx, istglob, ienglob,ap, bp, Tp, T, workarrp);

//    timestep_Euler(rank,size,nx,nxglob,istglob,ienglob,dt,dx,x,T,rhs,&t_comm);

   if((it%it_print == 0) && (it != 0))
    {
      if(rank==0) printf("Writing solution at time step no. %d, time = %lf\n", it, tcurr);
      get_exact_solution(nx,tcurr,x,Texactp);
      output_soln(rank,nx,it+1,tcurr,x,Tp,Texactp);
    }
  }
 MPI_Allgather(Texactp,nx,MPI_DOUBLE,Texact,nx,MPI_DOUBLE,MPI_COMM_WORLD);  
 MPI_Allgather(Tp,nx,MPI_DOUBLE,T,nx,MPI_DOUBLE,MPI_COMM_WORLD);
 err_norm(nxglob, nx, T, Texact, 2, &l2norm);
if(rank==0)
{
  printf("nx = %d; L2 norm (x-xex) = %.15e\n", nxglob, l2norm);
if(nxglob <= 10)
  {
    printf("\n*****Matrix A****\n");
//    print_matrix(nx, nx, a);
    printf("*****************\n\n");
    printf("*****Vector T****\n");
    print_vector(nxglob, T);
    printf("*****************\n\n");
    printf("*****Vector Tex****\n");
    print_vector(nxglob, Texact);
    printf("*****************\n\n");
    printf("*****Vector b****\n");
    //print_vector(nx, b);
    printf("*****************\n\n");
  }
}
 // output_error_norm(rank,size,nx,nxglob,tcurr,x,T,Texact);

//  free(rhs);
  free(T);
  free(x);
  free(Texact);

  MPI_Finalize();
  return 0;
}
