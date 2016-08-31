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
float distance_redshift(float z);

int main(int argc, char **argv)
{
  FILE *fp;
  char a,a1[10],a2[10],a3[10],aa[1000], fname[1000];
  int i, j, k,i1, j1, nmatched, idum, n1, n2, ngal;
  int **galid1, **galid2, xdum;
  float **galpos1, **galpos2, *redshift2, *mass2, *emass2; 
  float *magr, *magg, *ba90, *sersic_n;
  float *d4k, *d4k_err;
  int **ipos1, **ipos2, flag, nrepeat=0;
  float r, rmin, volume, vol_corr;
  float *mag_r, *mag_g, dx, dy, dz, x1, theta, phi, theta_max, rmax, *ra, *dec, *redshift, *mass, *rad, *sigma, *ang_rad;
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
  int cnt, n, imag;

  // for nonlocal
  float **group, *groupcnt, rr, rrsqr, rvir, xx[20];
  int *groupmem, igrp, nsample, i2, i3;

  if(argc>1)
    {
      imag = atoi(argv[1]);
      MAGNITUDE = -imag;
      MAGNITUDE2 = -imag;
      if(imag==18) REDSHIFT = 12000.0/SPEED_OF_LIGHT;
      if(imag==19) REDSHIFT = 19200.0/SPEED_OF_LIGHT;
      if(imag==20) REDSHIFT = 31000.0/SPEED_OF_LIGHT;
    }
  RADIUS = 10;
  if(argc>2)
    RADIUS = atof(argv[2]);
  fprintf(stderr,"Using RADIUS=%f\n",RADIUS);

  //read in all the galaxies from the VAGC
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/lss.dr72bright34.dat");
  ngal = filesize(fp);
  ra = vector(1,ngal);
  dec = vector(1,ngal);
  redshift = vector(1,ngal);
  indx = ivector(1,ngal);
  
  xg = vector(1,ngal);
  yg = vector(1,ngal);
  zg = vector(1,ngal);

  for(i=1;i<=ngal;++i) 
    {
      fscanf(fp,"%d %d %d %f %f %f",&idum,&idum,&idum,&ra[i],&dec[i],&redshift[i]);
      ra[i] *= PI/180.;
      dec[i] *= PI/180.;
      redshift[i] /= SPEED_OF_LIGHT;
      indx[i] = i;

      phi = ra[i];
      theta = PI/2.0 - dec[i];
      r = distance_redshift(redshift[i]);

      xg[i] = r*sin(theta)*cos(phi);
      yg[i] = r*sin(theta)*sin(phi);
      zg[i] = r*cos(theta);
		  
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in lss.dat\n");



  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/photo_evolve.dr72bright34.dat");
  //fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/bright1/photoinfo.dr72bright1.dat");
  if(filesize(fp)!=ngal)
    {
      fprintf(stderr,"ERROR: filesize mismatch for photinfo\n");
      exit(0);
    }
  mag_r = vector(1,ngal);
  mag_g = vector(1,ngal);
  for(i=1;i<=ngal;++i) 
    {
      fscanf(fp,"%d %f %f %f",&j,&x1,&mag_g[i],&mag_r[i]);
      //mag_r[i] = mag_r[i] + Q0*(1+Q1*(redshift[i]-QZ0))*(redshift[i]-QZ0);
      //mag_g[i] = mag_g[i] + Q0*(1+Q1*(redshift[i]-QZ0))*(redshift[i]-QZ0);
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in photoinfo.dat\n");


  // read in the halo mass of each galaxy
  group = matrix(1,ngal,1,4);
  for(i=1;i<=ngal;++i)
    group[i][1] = -1;

  groupmem = ivector(1,ngal);
  groupcnt = vector(1,ngal);

  for(i=1;i<=ngal;++i)
    groupmem[i] = groupcnt[i] = 0;

  sprintf(fname,"/Users/tinker/cosmo/CENSAT_DECOMP/DR7_RESULTS/clf_groups_M%d.prob",imag);
  fp = openfile(fname);
  nsample = filesize(fp);
  for(i=1;i<=nsample;++i)
    {
      fscanf(fp,"%s %d %d %d", aa, &i1,&i2,&i3);
      for(j=1;j<=6;++j)
	fscanf(fp,"%f",&xx[j]);
      fgets(aa,1000,fp);
      // i1 is the order in the vagc catalog
      groupmem[i1] = i3; // which group does this gal belong to
      // keep track of gal cnt within the defined sample.
      if(!(mag_r[i]>MAGNITUDE2 || redshift[i]>REDSHIFT || redshift[i]<CZMIN))
	groupcnt[i3]++;
      if(xx[2]>0.5)continue; // this gal is a satellite
      rvir = pow(3*xx[3]/(4*PI*2.775e11*0.25*200),1.0/3.0); //x[3] = mass
      group[i3][1] = rvir; // halo mass;
      group[i3][2] = xg[i1];
      group[i3][3] = yg[i1];
      group[i3][4] = zg[i1];
      //printf("%e %d %d %e %e\n",xx[3],i1,i3,ra[i1]*180/PI,dec[i1]*180/PI);
    }
  fclose(fp);

  volume =4./3.*PI*(pow(distance_redshift(REDSHIFT),3.0)-pow(distance_redshift(CZMIN),3.0))*
    9380.0/41253.0;//*0.59;//0.88;
  fprintf(stderr,"%e (Mpc/h)^3 --> [%e] per side\n",volume,pow(volume,0.3333));

  n=0;
  for(i=1;i<=ngal;++i)
    {
      if(mag_r[i]>MAGNITUDE2 || redshift[i]>REDSHIFT || redshift[i]<CZMIN)continue;
      n++;
    }
  mean_density = n/volume;
  RADIUS_SQR = RADIUS*RADIUS;
  mean_cnt = 4./3.*PI*RADIUS*RADIUS*RADIUS*mean_density;
  fprintf(stderr,"%d galaxies, mean in sphere= %f\n",n,mean_cnt);
  // exit(0);

  // find all the galaxies (in the sample) that are in the groups
  fprintf(stderr,"FInding group membership...\n");
  for(i=1;i<=ngal;++i)
    {
      if(group[i][1]<0)continue;
      cnt = 0;
      rr = group[i][1];
      rrsqr = rr*rr;
      for(j=1;j<=ngal;++j)
	{
	  if(mag_r[j]>MAGNITUDE)continue;
	  if(mag_r[j]>MAGNITUDE2)continue;
	  dx = fabs(group[i][2]-xg[j]);
	  if(dx>rr)continue;
	  dy = fabs(group[i][3]-yg[j]);
	  if(dy>rr)continue;
	  dz = fabs(group[i][4]-zg[j]);
	  if(dz>rr)continue;
	  
	  r = dx*dx + dy*dy + dz*dz;
	  if(r>rrsqr)continue;
	  cnt++;	  
	}
      groupcnt[i] = cnt;
    }
  fprintf(stderr,"Done FInding group membership...\n");

  for(i=1;i<=ngal;++i)
    {
      if(i%10000==0)fprintf(stderr,"%d\n",i);
      if(mag_r[i]>MAGNITUDE || redshift[i]>REDSHIFT || redshift[i]<CZMIN)continue;
      cnt = 0;
      igrp = groupmem[i];
      for(j=1;j<=ngal;++j)
	{
	  if(mag_r[j]>MAGNITUDE)continue;
	  if(mag_r[j]>MAGNITUDE2)continue;
	  //if(i==j)continue;
	  dx = fabs(group[igrp][2]-xg[j]);
	  if(dx>RADIUS)continue;
	  dy = fabs(group[igrp][3]-yg[j]);
	  if(dy>RADIUS)continue;
	  dz = fabs(group[igrp][4]-zg[j]);
	  if(dz>RADIUS)continue;

	  r = dx*dx + dy*dy + dz*dz;
	  if(r>RADIUS_SQR)continue;
	  cnt++;	  
	}
      // remove the number of galaxies within the group 
      cnt -= (groupcnt[igrp]+1);
      //density = cnt/mean_cnt;
      // remove volume of group as well
      density = cnt/(4./3.*PI*(RADIUS*RADIUS*RADIUS-pow(group[igrp][1],3.0)))/mean_density;
      printf("%e %d %f %f %f %d %.1f\n",density,cnt,mag_r[i],mag_g[i],redshift[i],igrp,groupcnt[igrp]);

    }
  
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

