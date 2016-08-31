#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

int main(int argc, char **argv)
{
  float **dbar, dx, dy, xmin, ymin, **xbar, **ybar, **ncnt;
  int nx,ny;
  float *dn4k, *Hdelta, *den;
  int i,j,n,ix,iy;
  FILE *fp;
  char aa[1000];

  xmin = -7;
  ymin = 0.8;
  dx = 0.5;
  dy = 0.05;
  nx = 30;
  ny = 30;

  ncnt = matrix(1,nx,1,ny);
  dbar = matrix(1,nx,1,ny);
  xbar = matrix(1,nx,1,ny);
  ybar = matrix(1,nx,1,ny);

  for(i=1;i<=nx;++i)
    for(j=1;j<=ny;++j)
      dbar[i][j] = xbar[i][j] = ybar[i][j] = ncnt[i][j] = 0;

  fp = openfile(argv[1]);
  n = filesize(fp)-12; //get rid of header
  
  dn4k = vector(1,n);
  Hdelta = vector(1,n);
  den = vector(1,n);
  
  fgets(aa,1000,fp); //header
  fgets(aa,1000,fp); //header

  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%f %f %f",&dn4k[i],&Hdelta[i],&den[i]);
      ix = (int)((Hdelta[i]-xmin)/dx)+1;
      iy = (int)((dn4k[i]-ymin)/dy)+1;
      //printf("%f %f %f %d %d\n",dn4k[i],Hdelta[i],den[i],ix,iy);
      if(ix<1 || ix>=nx)continue;
      if(iy<1 || iy>=ny)continue;
      dbar[ix][iy] += den[i];
      xbar[ix][iy] += Hdelta[i];
      ybar[ix][iy] += dn4k[i];
      ncnt[ix][iy]++;
    }

  for(i=1;i<=nx;++i)
    for(j=1;j<=ny;++j)
      {
	if(ncnt[i][j]>0)
	  printf("%e %e %e %.0f\n",xbar[i][j]/ncnt[i][j],ybar[i][j]/ncnt[i][j],dbar[i][j]/ncnt[i][j],ncnt[i][j]);
	else
	  printf("%f %ff -100 0\n",(i-0.5)*dx+xmin,(j-0.5)*dy+ymin);
      }


}
