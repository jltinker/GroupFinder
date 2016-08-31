#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define OMEGA_M 0.25
#define PI 3.14159
#define RHO_CRIT 2.775E+11
#define DELTA_HALO 200
#define SPEED_OF_LIGHT 3.0E+5
#define c_on_H0 3000.0
#define BIG_G 4.304E-9 /* BIG G in units of (km/s)^2*Mpc/M_sol */
#define G0 (1.0/sqrt(2.0*3.14159))
#define ROOT2 1.41421

//external functions
float qromo(float (*func)(float), float a, float b,
	     float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);


//local functions
float mass_from_luminosity(float x);
float angular_separation(float a1, float d1, float a2, float d2);
float distance_redshift(float z);

int main(int argc, char **argv)
{
  FILE *fp;
  char a,a1[10],a2[10],a3[10],aa[1000];
  int i, j, i1, j1, nmatched, idum, n1, n2, ngal;
  int **galid1, **galid2, xdum;
  float **galpos1, **galpos2, *redshift2, *mass2, *emass2; 
  float *magr, *magg, *ba90, *sersic_n;
  float *d4k, *d4k_err;
  int **ipos1, **ipos2, flag, nrepeat=0;
  float r, rmin;
  float *mag_r, *mag_g, dx, dy, dz, x1, theta, *ra, *dec, *redshift, *mass, *rad, *sigma, *ang_rad;
  int *icollided, *indx;

  float prob_rad, prob_ang, *prob_total;
  int imass;
  double nhalo[160]; 
  float nsat[160];
  double mbar[160];

  //testing vars
  float box_size = 384;
  float vx,vy,vz;

  for(i=0;i<160;++i)
    nsat[i] = nhalo[i] = mbar[i] = 0;

  //read in all the galaxies from the VAGC
  fp = openfile(argv[1]);
  ngal = filesize(fp);

  ra = vector(1,ngal);
  dec = vector(1,ngal);
  redshift = vector(1,ngal);
  sigma = vector(1,ngal);
  ang_rad = vector(1,ngal);
  indx = ivector(1,ngal);
  prob_total = vector(1,ngal);
  mass = vector(1,ngal);
  rad = vector(1,ngal);
  mag_r = vector(1,ngal);

  
  // read in galaxies from the mock
  // transform it into redshift space
  for(i=1;i<=ngal;++i) 
    {
      fscanf(fp,"%f %f %f %f %f %f",&ra[i],&dec[i],&redshift[i],&vx,&vy,&vz);
      redshift[i] += vz/100;
      if(redshift[i]>box_size)redshift[i] -= box_size;
      if(redshift[i]<=0)redshift[i] += box_size;
      indx[i] = i;

      redshift[i] *= 100.0; //convert to km/s

      fscanf(fp,"%f",&mag_r[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in lss.dat\n");

  // set up the mass/radius of each galaxy
  for(i=1;i<=ngal;++i)
    {
      mass[i] = mass_from_luminosity(mag_r[i]);
      rad[i] = pow(3*mass[i]/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
      //ang_rad[i] = rad[i]/distance_redshift(redshift[i]/SPEED_OF_LIGHT);
      ang_rad[i] = rad[i];
      sigma[i] = sqrt(BIG_G*mass[i]/rad[i]);
    }
  fprintf(stderr,"Done SHAMming the galaxies\n");

  //sort everything by luminosity
  sort2(ngal,mag_r,indx);
  fprintf(stderr,"Done sorting\n");

  
  for(i1=1;i1<=ngal;++i1)
    {
      i = indx[i1];
      if(mag_r[i1]>-20.999)
	prob_total[i] = -1;
      else
	prob_total[i] = 0;
    }
  //go through and find associated galaxies
  for(i1=1;i1<=ngal;++i1)
    {
      i = indx[i1];
      if(mag_r[i1]>-21.0)break;
      //printf("%f %f\n",mag_r[i1],redshift[i]/SPEED_OF_LIGHT);

      // only do things z<0.15 and M_r<-21
      if(i1%1000==0)fprintf(stderr,"%d\n",i1);
      if(mag_r[i1]>-21.0)continue;


      imass = (int)(log10(mass[i])*10.0);
      //mbar[imass] += mass[i];
      //nhalo[imass]++;
      // printf("%e %d\n",

      for(j1=i1+1;j1<=ngal;++j1)
	{
	  j = indx[j1];

	  // only do things z<0.15 and M_r<-21
	  if(mag_r[j1]>-21.0)continue;

	  //printf("%d %d %f %e %f %f\n",i1,j1,mag_r[i1],mass[i],rad[i],ang_rad[i]);

	  dx = fabs(ra[i]-ra[j]);
	  //printf("dx %e\n",dx);
	  if(dx>2*ang_rad[i])continue;
	  dy = fabs(dec[i]-dec[j]);
	  //printf("dy %e\n",dy);
	  if(dy>2*ang_rad[i])continue;
	  dz = fabs(redshift[i] - redshift[j]);
	  //printf("dz %e %e\n",dz, 6*sigma[i]);
	  if(dz>6*sigma[i])continue;

	  //theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);
	  theta = sqrt(dx*dx+dy*dy);
	  if(theta>ang_rad[i])continue;
	  //printf("HERE\n");

	  // check surface density of NFW halo
	  prob_ang = 1;

	  //bring in info on whether current galaxy is a satellite
	  prob_ang *= (1-prob_total[i]);
	  //if(prob_total[i]<0)printf("ERROR %f %d %d %f\n",prob_total[i],i,i1,mag_r[i1]);

	  // check radial probability
	  x1 = sqrt(sigma[i]*sigma[i] + 100.*100.*rad[i]*rad[i]);
	  //prob_rad = G0*exp(-dz*dz/(2*x1))/sqrt(x1);
	  prob_rad = erfc(dz/(ROOT2*x1));
	  //prob_rad = 0;
	  //if(dz<15*x1)prob_rad=1;


	  if(prob_total[j]<0)prob_total[j] = prob_ang*prob_rad;
	  else prob_total[j] += prob_ang*prob_rad;

	  nsat[imass]+=prob_ang*prob_rad;

	  //printf("%f %e %f %f %f %f %e %e %e\n",
	  //		 mag_r[i1],mass[i],rad[i],ang_rad[i],theta,mag_r[j1],mass[j],dz,x1);
	  //printf("%e\n",prob_total[j]);
	  //fflush(stdout);

	}


    }

  j=0;
  for(i=1;i<=ngal;++i)
    {
      if(prob_total[i]<0)continue;
      if(prob_total[i]>1)j++;
      if(prob_total[i]>1)prob_total[i]=1;
      imass = (int)(log10(mass[i])*10.0);
      mbar[imass] += mass[i]*(1-prob_total[i]);
      nhalo[imass] += (1-prob_total[i]);
    }
  fprintf(stderr,"%d\n",j);
  for(i=0;i<160;++i)
    if(nhalo[i]>0)
      printf("%e %e %e\n",mbar[i]/nhalo[i],1.0,nsat[i]/nhalo[i]);

}

// assuming that mass is in natural log.
float mass_from_luminosity(float x)
{
  FILE *fp;
  float y;
  int i;
  static float *mag, *mass, *zz;
  static int flag=1, n;

  if(flag)
    {
      flag = 0;
      fp = openfile("sham_lum2mass.dat");
      n = filesize(fp);
      mag = vector(1,n);
      mass = vector(1,n);
      zz = vector(1,n);
      for(i=1;i<=n;++i)
	fscanf(fp,"%f %f",&mag[i],&mass[i]);
      fclose(fp);
      spline(mag,mass,n,1.0E+30,1.0E+30,zz);
      fprintf(stderr,"Done setting up mag2mass\n");
    }

  splint(mag,mass,zz,n,x,&y);  
  return exp(y);
}


float func_dr1(float z)
{
  return pow(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M),-0.5);
}
float distance_redshift(float z)
{
  float x;
  if(z<=0)return 0;
  //fprintf(stderr,"blah %f\n",z);
  x= c_on_H0*qromo(func_dr1,0.0,z,midpnt);
  //fprintf(stderr,"%f %f\n",z,x);
  return x;
}

float angular_separation(float a1, float d1, float a2, float d2)
{
  float cd1,cd2,sd1,sd2,ca1a2,sa1a2;

  return atan((sqrt(cos(d2)*cos(d2)*sin(a2-a1)*sin(a2-a1) + 
		    pow(cos(d1)*sin(d2) - sin(d1)*cos(d2)*cos(a2-a1),2.0)))/
	      (sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a2-a1)));
}
