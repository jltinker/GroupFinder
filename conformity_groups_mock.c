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
  float *mag_r, *mag_g, dx, dy, dz, x1, theta, phi, theta_max, rmax, *ra, *dec, *redshift, *mass, *rad, *sigma, *ang_rad, *mgal, *grcolor, *colorgrp;
  int *icollided, *indx, galid;

  float prob_rad, prob_ang, *prob_total, p, meanmass;
  int imass, ncnt2; 
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
    *fquench,*fqsat,*fqcen,*dn4kgrp;
  float ANGLE, DELTAZ, MAGLO, MAGHI, MLO, PLIM, DELTA_M, BREAK;
  float fq_jack[50][50], fqs_jack[50][50], fqc_jack[50][50], galcnt_jack[50][50], satcnt_jack[50][50];

  int BLUE;
  float CENTHRESH;

  int nby5, nsby5, *subindx, *ijack, ii, id;
  float err_qc, err_qs, err_qt;
  char fname[1000];
  
  for(i=0;i<50;++i)
    for(j=0;j<50;++j)
      fq_jack[i][j] = fqs_jack[i][j] = fqc_jack[i][j] = galcnt_jack[i][j] = satcnt_jack[i][j] = 0;

  DELTAZ = 1000;
  ANGLE = 10*PI/180.;
  MAGLO = -19.0;
  MAGHI = -20.0;
  MLO = 9.6;
  DELTA_M = 0.2;
  PLIM = 0.5;


  if(argc>1)
    {
      imag = atoi(argv[1]);
    }
  MLO = 9.6;
  if(argc>2)
    MLO = atof(argv[2]);
  CENTHRESH = 0.5;
  if(argc>3)
    CENTHRESH = atof(argv[3]);
  BLUE = 0;
  if(argc>4)
    BLUE = atoi(argv[4]);
  DELTA_M = 0.2;
  if(argc>5)
    DELTA_M = atof(argv[5]);
  if(argc>6)
    BREAK= atof(argv[6]);

  fprintf(stderr,"input: %d %f %f %d %f\n",imag, MLO, CENTHRESH, BLUE, DELTA_M);

  sprintf(fname,"clf_groups_M19_%d.galdata",imag);
  fp = openfile(fname);
  ngal = filesize(fp);
  ra = vector(1,ngal);
  dec = vector(1,ngal);
  redshift = vector(1,ngal);
  mgal = vector(1,ngal);
  indx = ivector(1,ngal);
  subindx = ivector(1,ngal);
  ijack = ivector(1,ngal);

  xg = vector(1,ngal);
  yg = vector(1,ngal);
  zg = vector(1,ngal);

  mag_r = vector(1,ngal);
  mag_g = vector(1,ngal);
  dn4k = vector(1,ngal);
  grcolor = vector(1,ngal);

  for(i=1;i<=ngal;++i) 
    {
      fscanf(fp,"%s %d %f %f %f %f %f %f",aa,&j,&ra[i],&dec[i],&redshift[i],&x1,&mag_r[i],&dn4k[i]);
      mgal[i] = fabs(mag_r[i]);
      indx[i] = i;

      phi = ra[i];
      theta = PI/2.0 - dec[i];
      r = redshift2distance(redshift[i]);

      xg[i] = r*sin(theta)*cos(phi);
      yg[i] = r*sin(theta)*sin(phi);
      zg[i] = r*cos(theta);

      xg[i] = dec[i];

      fgets(aa,1000,fp);
      //printf("BOO %f %f %f\n",mag_r[i],mgal[i],dn4k[i]);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in lss.dat\n");


  // lets get a jackknife sampling index
  fprintf(stderr,"jackknifing sample...\n");
  
  // sort everything in dec (xg is not actually used)
  sort2(ngal, xg, indx);

  // now go in quintiles of dec and sort by ra.
  nby5 = ngal/5;
  for(i=0;i<5;++i)
    {
      for(n=0,j=i*nby5+1;j<=(i+1)*nby5;++j)
	{
	  n++;
	  id = indx[j];
	  yg[n] = ra[id];
	  subindx[n] = id;
	}
      // sort by ra
      sort2(n,yg,subindx);
      // now divide these into groups of five
      nsby5 = n/5;
      for(ii=0;ii<5;++ii)
	for(k=0,j=ii*nsby5+1;j<=(ii+1)*nsby5;++j)
	  {
	    k++;
	    id = subindx[j];
	    ijack[id] = i + ii*5;
	  }
    }

  for(i=1;i<=-ngal;++i)
    printf("BOO %d %f %f\n",ijack[i],ra[i],dec[i]);

  //read in the probabilities
  psat = vector(1,ngal);
  sprintf(fname,"clf_groups_M19_%d.prob",imag);
  fp = openfile(fname);

  n = filesize(fp);
  if(n!=ngal)
    {
      fprintf(stderr,"ERROR: filesize mismatch.\n");
      exit(0);
    }
  for(i=1;i<=n;++i)
    {
      fscanf(fp,"%s %d %d %d %f %f",aa,&i1,&j,&i1,&x1,&x2);
      fgets(aa,1000,fp);
      psat[i] = x2;
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
  /*
  if(imag==1) 
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M19_M9.8.groups");
  if(imag==2)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18.5_M9.6.groups");
  if(imag==3)
    fp = openfile("/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M18_M9.4.groups");
  ngrp = filesize(fp);
  */

  // don't do the groups-- just do the central galaxies
  ngrp = 0;
  for(i=1;i<=ngal;++i)
    if(psat[i]<CENTHRESH)ngrp++;
  fprintf(stderr,"numgrps: %d\n",ngrp);

  mgrp = vector(1,ngrp);
  ragrp = vector(1,ngrp);
  decgrp = vector(1,ngrp);
  zgrp = vector(1,ngrp);
  rcogrp = vector(1,ngrp);
  radgrp = vector(1,ngrp);
  dn4kgrp = vector(1,ngrp);
  colorgrp = vector(1,ngrp);

  // temp
  for(i=1;i<=ngal;++i)
    xg[i] = ijack[i];

  i=0;
  for(j=1;j<=ngal;++j)
    {
      if(psat[j]>=CENTHRESH)continue;
      i++;
      mgrp[i] = (mgal[j]); // now is the central gal mass
      ragrp[i] = ra[j];
      decgrp[i] = dec[j];
      zgrp[i] = redshift[j];
      rcogrp[i] = redshift2distance(zgrp[i]);
      radgrp[i] = pow(1.0E12*3/(4*PI*200*2.775e11*OMEGA_M),1.0/3.0); //whatever
      dn4kgrp[i] = dn4k[j];
      colorgrp[i] = grcolor[j];
      ijack[i] = xg[j];
    }
  


  // log bins
  nbins = 17;
  rmin = 0.1;
  dlogr = log(10/0.1)/nbins;

  //linear bins
  nbins = 10;
  rmin = 0.1;
  dlogr = 1;
  galcnt = ivector(1,nbins);
  rbar = vector(1,nbins);
  satcnt = ivector(1,nbins);
  fquench = vector(1,nbins);
  fqsat = vector(1,nbins);
  fqcen = vector(1,nbins);

  for(i=1;i<=nbins;++i)
    fqsat[i] = fqcen[i] = fquench[i] = galcnt[i] = rbar[i] = satcnt[i] = mbar[i] = 0;
  
  ncnt = 0;
  ncnt2 = 0;
  meanmass = 0;
  for(j=1;j<=ngrp;++j)
    {
      // only do centrals, but make it fixed stellar mass, not halo mass
      if(mgrp[j]<MLO || mgrp[j]>MLO+DELTA_M)continue;
      // decision tree on central Dn400
      if(BLUE) {
	if(dn4kgrp[j]>BREAK)continue; }
      else {
	if(dn4kgrp[j]<BREAK)continue; }
      ncnt++;

      // also only choose galaxies in the same mass bin
      for(i=1;i<=ngal;++i)
	{	  
	  if(mgal[i]<MLO || mgal[i]>MLO+DELTA_M)continue;

	  dz = fabs(redshift[i]-zgrp[j])*SPEED_OF_LIGHT;
	  if(dz>DELTAZ)continue;
	  dx = fabs(ra[i]-ragrp[j]);
	  if(dx>ANGLE)continue;
	  dy = fabs(dec[i]-decgrp[j])*cos(dec[i]); 
	  if(dy>ANGLE)continue;
	  theta = angular_separation(ra[i],dec[i],ragrp[j],decgrp[j]);
	  if(theta>ANGLE)continue;

	  dr = theta*rcogrp[j];
	  //ibin = log(dr/rmin)/dlogr + 1;
	  if(dr<rmin)continue;
	  ibin = dr/dlogr + 1;

	  if(ibin>nbins)continue;
	  if(ibin<=0)continue;

	  galcnt[ibin]++;
	  rbar[ibin]+=dr;
	  if(psat[i]>PLIM)
	    satcnt[ibin]++;
	  
	  if(psat[i]<PLIM)
	    mbar[ibin] += pow(10.0,mgal[i]);
	  if(dn4k[i]>BREAK)
	    fquench[ibin]++;
	  if(dn4k[i]>BREAK && psat[i]>PLIM)
	    fqsat[ibin]++;
	  if(dn4k[i]>BREAK && psat[i]<=PLIM)
	    fqcen[ibin]++;

	  // jackknife sample on the primary
	  ii = ijack[j];
	  galcnt_jack[ibin][ii]++;
	  if(psat[i]>PLIM)
	    satcnt_jack[ibin][ii]++;
	  if(dn4k[i]>BREAK)
	    fq_jack[ibin][ii]++;
	  if(dn4k[i]>BREAK && psat[i]>PLIM)
	    fqs_jack[ibin][ii]++;
	  if(dn4k[i]>BREAK && psat[i]<=PLIM)
	    fqc_jack[ibin][ii]++;


	}
    }
  fprintf(stderr,"%d clusters\n",ncnt);

  for(i=1;i<=nbins;++i)
    if(galcnt[i])
      {
	ncnt = galcnt[i]-satcnt[i];

	// loop over the jacks
	err_qt = err_qs = err_qc = 0;
	for(j=0;j<25;++j)
	  {
	    x1 = fquench[i]/galcnt[i];
	    x2 = (fquench[i]-fq_jack[i][j])/(galcnt[i]-galcnt_jack[i][j]);
	    err_qt += (x2-x1)*(x2-x1);

	    x1 = fqsat[i]/satcnt[i];
	    x2 = (fqsat[i]-fqs_jack[i][j])/(satcnt[i]-satcnt_jack[i][j]);
	    err_qs += (x2-x1)*(x2-x1);

	    x1 = fqcen[i]/(galcnt[i]-satcnt[i]);
	    x2 = (fqcen[i]-fqc_jack[i][j])/((galcnt[i]-satcnt[i]) - (galcnt_jack[i][j]-satcnt_jack[i][j]));
	    err_qc += (x2-x1)*(x2-x1);
	  }
	err_qt = sqrt(25./24.*err_qt);
	err_qc = sqrt(25./24.*err_qc);
	err_qs = sqrt(25./24.*err_qs);

	if(galcnt[i]>satcnt[i])
	  {
	    printf("%d %f %f %f %f %d %d %f %f %f %f\n",i,rbar[i]/galcnt[i],fquench[i]/galcnt[i],
		   fqsat[i]/satcnt[i],fqcen[i]/(galcnt[i]-satcnt[i]),satcnt[i],
		   galcnt[i],log10(mbar[i]/ncnt),err_qt, err_qs, err_qc);
	  }
	else
	  {
	    printf("%d %f %f %f %f %d %d %f\n",i,rbar[i]/galcnt[i],fquench[i]/galcnt[i],
		   fqsat[i]/satcnt[i],0.0,satcnt[i],galcnt[i],log10(mbar[i]/ncnt));
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

