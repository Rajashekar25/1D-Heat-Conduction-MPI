#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void conjug_grad(int nx, int nxloc, int istglob, int ienglob, double **ap, double *bp, double *xp, double *x, double *workarrp)
{
  int i, j, iter, max_iter;
  double *rkp, *rk, *pkp, *pk, *vec2p, *vec2, locnumer_alp, numer_alp, locnumer_bet, numer_bet, locdenom, denom, alpk, betk, locnormb, normb, tol;
  max_iter = 10*nx; tol = 1.0e-12;
 // x = (double *)malloc(nx*sizeof(double));
  rk = (double *)malloc(nx*sizeof(double));
  pk = (double *)malloc(nx*sizeof(double));
//  vec2 = (double *)malloc(nx*sizeof(double));
  vec_norm(nx, nxloc, bp, &locnormb);
 MPI_Allreduce(&locnormb, &normb,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  normb = sqrt(normb/(double)nx);
  rkp = &workarrp[0]; pkp = &workarrp[nxloc]; vec2p = &workarrp[2*nxloc];
  
  // initialization
 MPI_Allgather(xp,nxloc,MPI_DOUBLE,x,nxloc,MPI_DOUBLE,MPI_COMM_WORLD);
  matrix_multiply(nx, nxloc, istglob, ienglob, ap, xp, rkp, x);
  for(i=0; i<nxloc; i++)
  {
    rkp[i] = bp[i] - rkp[i];
    pkp[i] = rkp[i];
  }
  MPI_Allgather(pkp,nxloc,MPI_DOUBLE,pk,nxloc,MPI_DOUBLE,MPI_COMM_WORLD);
  matrix_multiply(nx, nxloc, istglob, ienglob, ap, pkp, vec2p, pk);

  dot_prod(nxloc, istglob, ienglob, rkp, rkp, &locnumer_alp);
  MPI_Allreduce(&locnumer_alp, &numer_alp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(iter=0; iter<max_iter; iter++)
  {
    dot_prod(nxloc, istglob, ienglob, pkp, vec2p, &locdenom);
   MPI_Allreduce(&locdenom, &denom, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(fabs(denom)>1.0e-24)
      alpk = numer_alp/denom;
    else
      printf("\n Iter no. %d; Residual %.7e is very close to zero. Check details.", iter, denom);

   for(i=0; i<nxloc; i++)
   {
     xp[i] += alpk*pkp[i];
     rkp[i] -= alpk*vec2p[i];
   }

   dot_prod(nxloc, istglob, ienglob, rkp, rkp, &locnumer_bet);
  MPI_Allreduce(&locnumer_bet, &numer_bet, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   betk = numer_bet/numer_alp;

   for(i=0; i<nxloc; i++)
     pkp[i] = rkp[i] + betk*pkp[i];

   // stopping condition
   if(fabs(sqrt(numer_bet/(double)nx)/normb) < tol)
     break;

   // prepare for the next iteration
   numer_alp = numer_bet;
   MPI_Allgather(pkp,nxloc,MPI_DOUBLE,pk,nxloc,MPI_DOUBLE,MPI_COMM_WORLD);
   matrix_multiply(nx, nxloc, istglob, ienglob, ap, pkp, vec2p, pk);
  }
}