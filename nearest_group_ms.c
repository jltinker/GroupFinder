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
#define Q0 2.0
#define Q1 -1.0
#define QZ0 0.1
#define THIRD (1.0/3.0)

//numerical recipes
float qromo(float (*func)(float), float a, float b,
	     float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);

//external functions
float density2halo(float galaxy_density);
float halo_abundance(float m);
float subhalo_abundance(float m);
float func_nhalo(float m);

//local functions
float mass_from_luminosity(float x);
float angular_separation(float a1, float d1, float a2, float d2);
float redshift2distance(float z);

int main(int argc, char **argv)
{
  FILE *fp;
  char a,a1[10],a2[10],a3[10],aa[1000];
  int i, j, k,i1, j1, nmatched, idum, n1, n2, ngal;
  int **galid1, **galid2, xdum;
  float **galpos1, **galpos2, *redshift2, *mass2, *emass2; 
  float *magr, *magg, *ba90, *sersic_n;
  float *d4k, *d4k_err;
  int **ipos1, **ipos2, flag, nrepeat=0;
  float r, volume, vol_corr;
  float *mag_r, *mag_g, dx, dy, dz, x1, theta, phi, theta_max, rmax, *ra, *dec, *redshift, *mass, *rad, *sigma, *ang_rad, *mgal;
  int *icollided, *indx;

  float prob_rad, prob_ang, *prob_total, p;
  int imass; 
  double nsat[160], nhalo[160], nsat_cur, nsub[160]; 
  double nsat2[160], nhalo2[160], mbar2[160]; 
  double mbar[160];

  float REDSHIFT=19200.0/SPEED_OF_LIGHT, MAGNITUDE=-19.0, CZMIN=7000.0/SPEED_OF_LIGHT, CZBUF=0;
  //float REDSHIFT=31000.0/SPEED_OF_LIGHT, MAGNITUDE=-20.0, CZMIN=7000.0/SPEED_OF_LIGHT, CZBUF=0;
  //float REDSHIFT=12000.0/SPEED_OF_LIGHT, MAGNITUDE=-18.0, CZMIN=7000.0/SPEED_OF_LIGHT, CZBUF=0;
  float MAGNITUDE2 = -19;

  // for environoments
  float *xg,*yg,*zg;
  float mean_density, mean_cnt;
  float density;
  float RADIUS, RADIUS_SQR;
  int cnt, n, imag, ncnt;

  // for nearest group
  float *rbar, dlogr, rmin;
  int *galcnt, *satcnt, nbins, ngrp, ibin;
  float *psat, *rcogrp, *mgrp, *ragrp, *decgrp, *zgrp, xx[100], x2, dr, *radgrp, *dn4k,
    *fquench,*fqsat,*fqcen;
  float ANGLE, DELTAZ, MAGLO, MAGHI, MLO, PLIM, DELTA_M;

  float *groupsep, sigmav, *mhost;

  DELTAZ = 2000;
  ANGLE = 10*PI/180.;
  MAGLO = -19.0;
  MAGHI = -20.0;
  MLO = 14.0;
  DELTA_M = 2.0;
  PLIM = 0.5;


  if(argc>1)
    {
      imag = atoi(argv[1]);
    }
  if(imag==1)
    {
      MAGLO = -9.8;
      MAGHI = -10.2;
    }
  if(imag==2)
    {
      MAGLO = -9.6;
      MAGHI = -9.8;
    }
  if(imag==3)
    {
      MAGLO = -9.4;
      MAGHI = -9.8;
    }
  if(imag==0)
    {
      imag = 1;
      MAGLO = -10.0;
      MAGHI = -10.4;
    }
  if(imag==-1)
    {
      imag = 1;
      MAGLO = -10.2;
      MAGHI = -10.4;
    }



  //read in all the galaxies from the VAGC
  if(imag==0)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M19_M9.8.galdata_corr");
  if(imag==1)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M19_M9.8.galdata_corr");
  if(imag==2)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18.5_M9.6.galdata_corr");
  if(imag==3)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18_M9.4.galdata_corr");
  ngal = filesize(fp);
  ra = vector(1,ngal);
  dec = vector(1,ngal);
  redshift = vector(1,ngal);
  mgal = vector(1,ngal);
  indx = ivector(1,ngal);
  
  xg = vector(1,ngal);
  yg = vector(1,ngal);
  zg = vector(1,ngal);

  mag_r = vector(1,ngal);
  mag_g = vector(1,ngal);
  dn4k = vector(1,ngal);
  mhost = vector(1,ngal);

  for(i=1;i<=ngal;++i) 
    {
      fscanf(fp,"%s %d %f %f %f %f %f %f %f %f %f",aa,&idum,&mag_r[i],&mag_g[i],&redshift[i],&dn4k[i],
	     &x1,&x1,&mgal[i],&ra[i],&dec[i]);      
      indx[i] = i;

      phi = ra[i];
      theta = PI/2.0 - dec[i];
      r = redshift2distance(redshift[i]/SPEED_OF_LIGHT);

      xg[i] = r*sin(theta)*cos(phi);
      yg[i] = r*sin(theta)*sin(phi);
      zg[i] = r*cos(theta);

      mag_r[i] = -log10(mgal[i]);
		  

      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in lss.dat\n");


  //read in the probabilities
  psat = vector(1,ngal);
  fprintf(stderr,"here\n");
  if(imag==1) 
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M19_M9.8.prob");
  if(imag==2)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18.5_M9.6.prob");
  if(imag==3)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18_M9.4.prob");
  fprintf(stderr,"here\n");
  n = filesize(fp);
  fprintf(stderr,"here\n");
  if(n!=ngal)
    {
      fprintf(stderr,"ERROR: filesize mismatch.\n");
      exit(0);
    }
  fprintf(stderr,"here\n");
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%s %d %d %d",aa,&i1,&j,&i1);
      for(j=1;j<=9;++j)
	fscanf(fp,"%f",&xx[j]);
      psat[i] = xx[2];
      mhost[i] = log10(xx[9]);
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done with probabilites\n");

  n=0;
  for(i=1;i<=ngal;++i)
    {
      if(mag_r[i]>MAGNITUDE2 || redshift[i]>REDSHIFT || redshift[i]<CZMIN)continue;
      n++;
    }

  // read in all the groups
  if(imag==1) 
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M19_M9.8.groups");
  if(imag==2)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18.5_M9.6.groups");
  if(imag==3)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18_M9.4.groups");
  //fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M19.groups");
  ngrp = filesize(fp);
  mgrp = vector(1,ngrp);
  ragrp = vector(1,ngrp);
  decgrp = vector(1,ngrp);
  zgrp = vector(1,ngrp);
  rcogrp = vector(1,ngrp);
  radgrp = vector(1,ngrp);

  for(i=1;i<=ngrp;++i)
    {
      fscanf(fp,"%s %d %d",aa,&j,&j);
      for(j=1;j<=9;++j)
	fscanf(fp,"%f",&xx[j]);
      mgrp[i] = log10(xx[1]);
      ragrp[i] = xx[7];
      decgrp[i] = xx[8];
      zgrp[i] = xx[9];
      rcogrp[i] = redshift2distance(zgrp[i]/SPEED_OF_LIGHT);
      radgrp[i] = pow(xx[1]*3/(4*PI*200*2.775e11*OMEGA_M),1.0/3.0);
    }
  fclose(fp);
  fprintf(stderr,"Done reading [%d] lines from groups files\n",ngrp);
  

  nbins = 17;
  rmin = 0.1;
  dlogr = log(10/0.1)/nbins;
  galcnt = ivector(1,nbins);
  rbar = vector(1,nbins);
  satcnt = ivector(1,nbins);
  fquench = vector(1,nbins);
  fqsat = vector(1,nbins);
  fqcen = vector(1,nbins);

  for(i=1;i<=nbins;++i)
    fqsat[i] = fqcen[i] = fquench[i] = galcnt[i] = rbar[i] = satcnt[i] = 0;
  
  // vector for all galaxies
  groupsep = vector(1,ngal);
  for(i=1;i<=ngal;++i)
    groupsep[i] = 10;

  ncnt = 0;
  for(j=1;j<=ngrp;++j)
    {
      //if(mgrp[j]<MLO || mgrp[j]>MLO+DELTA_M)continue;
      ncnt++;

      sigmav = sqrt(BIG_G*pow(10.0,mgrp[j])/2.0/radgrp[j]*(1+redshift[j]/SPEED_OF_LIGHT));
      
      for(i=1;i<=ngal;++i)
	{	  
	  //if(mag_r[i]>MAGLO || mag_r[i]<MAGHI)continue;

	  dz = fabs(redshift[i]-zgrp[j]);
	  if(dz>DELTAZ)continue;
	  dx = fabs(ra[i]-ragrp[j]);
	  if(dx>ANGLE)continue;
	  dy = fabs(dec[i]-decgrp[j]);
	  if(dy>ANGLE)continue;
	  theta = angular_separation(ra[i],dec[i],ragrp[j],decgrp[j]);
	  if(theta>ANGLE)continue;

	  

	  dr = theta*rcogrp[j]/radgrp[j];
	  if(dr<groupsep[i] && dz<2.5*sigmav && psat[i]<0.5 && dr>1 && mhost[i]<mgrp[j]) { 
	    // printf("%d %d %e %e %e %e %e %e %e %e %e\n",j,i,
	    //	   dr,theta,rcogrp[j],radgrp[j],mgrp[j],dz,sigmav,redshift[i],zgrp[j]);
	    groupsep[i] = dr; }
	  ibin = log(dr/rmin)/dlogr + 1;

	  if(ibin>nbins)continue;
	  if(ibin<=0)continue;

	  galcnt[ibin]++;
	  rbar[ibin]+=dr;
	  if(psat[i]>PLIM)
	    satcnt[ibin]++;
	  if(dn4k[i]>1.6)
	    fquench[ibin]++;
	  if(dn4k[i]>1.6 && psat[i]>PLIM)
	    fqsat[ibin]++;
	  if(dn4k[i]>1.6 && psat[i]<=PLIM)
	    fqcen[ibin]++;
	}
    }
  fprintf(stderr,"%d clusters\n",ncnt);

  for(i=1;i<=ngal;++i)
    printf("%f %f %f\n",groupsep[i],psat[i],fabs(mag_r[i]));
  exit(0);

  for(i=1;i<=nbins;++i)
    if(galcnt[i])
      {
	if(galcnt[i]>satcnt[i])
	  {
	    printf("%d %f %f %f %f %d %d\n",i,rbar[i]/galcnt[i],fquench[i]/galcnt[i],
		   fqsat[i]/satcnt[i],fqcen[i]/(galcnt[i]-satcnt[i]),satcnt[i],galcnt[i]);
	  }
	else
	  {
	    printf("%d %f %f %f %f %d %d\n",i,rbar[i]/galcnt[i],fquench[i]/galcnt[i],
		   fqsat[i]/satcnt[i],0.0,satcnt[i],galcnt[i]);
	  }
      }	  
}

float func_dr1(float z)
{
  return pow(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M),-0.5);
}
float redshift2distance(float z)
{
  float x;
  if(z<=0)return 0;
  x= c_on_H0*qromo(func_dr1,0.0,z,midpnt);
  return x;
}

float angular_separation(float a1, float d1, float a2, float d2)
{
  float cd1,cd2,sd1,sd2,ca1a2,sa1a2;

  return atan((sqrt(cos(d2)*cos(d2)*sin(a2-a1)*sin(a2-a1) + 
		    pow(cos(d1)*sin(d2) - sin(d1)*cos(d2)*cos(a2-a1),2.0)))/
	      (sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a2-a1)));
}

