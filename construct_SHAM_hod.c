
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define csign(x) (x < 0.0 ? -1 : 1)
#define ind(a,b,c) (a)*n*n+(b)*n+(c)
#define indfft(a,b,c) (a)*(n+2)*n+(b)*(n+2)+(c)
#define PI 3.1415926535898

void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart);


int main(int argc, char *argv[])
{
  float rmin,rmax,xv,yv,zv,rcube,dx,dy,dz,rhalf,r,rp,x1,dlogr,delta,
    delta0=-0.8,meanmass[100],delta1=0.8,meanmass1,radii[100],radii2[100],vv,pmass,rad;
  int i,j,k,nr,np,nv,iv[100],ilo=1920,io[100],*imass;
  double iu[100];

  FILE *fp,*fp2;
  float *rsqr,*xg,*yg,*zg,*vxg,*vyg,*vzg,*x,*y,*z,*fdat,znow,x2,*mass,*rvir,*hostmass;
  int *meshparts, ***meshstart,nmesh,meshfac,nbrmax,*indx,ngal,nlim,nrad,*idat;
  char aa[1000];
  int *isat,*haloid;
  

  float *den,rscale,redshift,h0;
  int n,n3p,ldx,ip,jp,kp,i1,j1,k1,dark_matter=0,nmax;

  if(argc<3)
    {
      fprintf(stderr," ./vpf halofile boxsize > outfile\n\n");
      exit(0);
    }

  rcube = atof(argv[2]);
  rhalf=rcube*0.5;

  if(argc>4)
    dark_matter=atoi(argv[4]);

  /* Read in the galaxies from the galaxy file.
   */ 
  fp=openfile(argv[1]);
  ngal=filesize(fp)-5; //5 is the buffer

  xg=malloc(ngal*sizeof(float));
  yg=malloc(ngal*sizeof(float));
  zg=malloc(ngal*sizeof(float));
  /*
  vxg=malloc(ngal*sizeof(float));
  vyg=malloc(ngal*sizeof(float));
  vzg=malloc(ngal*sizeof(float));
  */
  mass=malloc(ngal*sizeof(float));
  hostmass=malloc(ngal*sizeof(float));
  rvir=malloc(ngal*sizeof(float));
  haloid = malloc(ngal*sizeof(int));
  isat=malloc(ngal*sizeof(int));
  rewind(fp);

  //get rid of the header
  for(i=1;i<=5;++i)
    fgets(aa,1000,fp);

  for(x1=k=i=0;i<ngal;++i)
    {
      fscanf(fp,"%f %f %f %f %f %f %f",&xg[i],&yg[i],&zg[i],&x1,&x1,&x1,&mass[i]);
      mass[i] *= 1.3;
      rvir[i] = pow(3*mass[i]/(4*PI*2.775e11*200*0.25),1.0/3.0);
      haloid[i] = i;
      isat[i] = -1;
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Read [%d] galaxies from [%s]\n",ngal,argv[1]);

  //sort2(ngal,&mass[0],&haloid[0]);

  indx=malloc(ngal*sizeof(int));
  rsqr=malloc(ngal*sizeof(float));
  fprintf(stderr,"starting meshlink for galaxies.\n");
  nmesh=0;
  rad = rvir[0]*1.3;
  meshlink2(ngal,&nmesh,0.0,rcube,rad,xg,yg,zg,&meshparts,&meshstart,meshfac);
  fprintf(stderr,"done with meshlink\n");
  
  for(i1=0;i1<ngal;++i1)
    {
      // if(i1%10000==0)fprintf(stderr,"%d\n",i1);
      i = haloid[i1];
      if(isat[i]>=0)continue;

      nbrmax=ngal;
      nbrsfind2(0.0,rcube,rvir[i],nmesh,xg[i],yg[i],zg[i],&nbrmax,indx,rsqr,xg,yg,zg,
		meshparts,meshstart);

      for(j=0;j<nbrmax;++j)
	{
	  if(indx[j]==i)continue;
	  if(isat[indx[j]]!=-1)continue;
	  isat[indx[j]]=i;
	  hostmass[indx[j]]=mass[i];
	}
    }
  for(i1=0;i1<ngal;++i1)
    {
      i = haloid[i1];
      if(isat[i]<0)hostmass[i]=mass[i1];
      printf("%e %d %e\n",mass[i1],isat[i],hostmass[i]);
    }


}

