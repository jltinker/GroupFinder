#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"

int main(int argc, char **argv)
{
  float *x1, *x3, *x4, *x5;
  int *imagr, *ims, i, j, k, nmagr, nms, *ix2;
  char fname[1000], aa[1000];
  FILE *fpr, *fps, *fpa;

  fpr = openfile(argv[1]);
  fpa = openfile(argv[2]);
  nmagr = filesize(fpr);
  imagr = ivector(1,nmagr);
  ix2 =  ivector(1,nmagr);
  x1 = vector(1,nmagr);
  x3 = vector(1,nmagr);
  x4 = vector(1,nmagr);
  x5 = vector(1,nmagr);

 
  for(i=1;i<=nmagr;++i)
    {
      fscanf(fpr,"%s %d",aa,&imagr[i]);
      fgets(aa,1000,fpr);
      fscanf(fpa,"%f %d %f %f %f",&x1[i],&ix2[i],&x3[i],&x4[i],&x5[i]);
    }

  fps = openfile(argv[3]);
  nms = filesize(fps);
  ims = ivector(1,nms);

  for(i=1;i<=nms;++i)
    {
      fscanf(fps,"%s %d",aa,&ims[i]);
      fgets(aa,1000,fps);
    }
  fclose(fpr);
  fclose(fps);
  fclose(fpa);
  
  for(i=1;i<=nms;++i)
    {
      for(j=1;j<=nmagr;++j)
	if(imagr[j] == ims[i])break;
      if(j>nmagr)
	{
	  fprintf(stderr,"ERROR: no match for %d\n",i);
	  exit(0);
	}
      printf("%e %d %f %f %f\n",x1[j],ix2[j],x3[j],x4[j],x5[j]);
    }


}
