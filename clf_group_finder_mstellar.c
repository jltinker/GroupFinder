#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define OMEGA_M 0.25
#define PI 3.141592741
#define RHO_CRIT 2.775E+11
#define DELTA_HALO 200
#define SPEED_OF_LIGHT 3.0E+5
#define c_on_H0 2997.92
#define BIG_G 4.304E-9 /* BIG G in units of (km/s)^2*Mpc/M_sol */
#define G0 (1.0/sqrt(2.0*3.14159))
#define ROOT2 1.41421
#define Q0 2.0
#define Q1 -1.0
#define QZ0 0.1
#define THIRD (1.0/3.0)
#define ANG (PI/180.0)
#define RT2PI 2.50663

#define CZMIN 7000.0
#define CZBUF 0
//#define REDSHIFT (12000.0/SPEED_OF_LIGHT)// -18
//#define MAGNITUDE -18.0
//#define REDSHIFT (19200.0/SPEED_OF_LIGHT)// -19
//#define MAGNITUDE -19.0
//#define REDSHIFT (31000.0/SPEED_OF_LIGHT)// -20
//#define MAGNITUDE -20.0

//numerical recipes
float qromo(float (*func)(float), float a, float b,
	     float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void sort2(int n, float arr[], int id[]);
float qtrap(float (*func)(float), float a, float b);

//external functions
float density2halo(float galaxy_density);
float density2host_halo(float galaxy_density);
float halo_abundance(float m);
float subhalo_abundance(float m);
float subhalo_abundance2(float m);
float func_nhalo(float m);
float func_subhalo1(float mhost);
void scatter_test(void);
float most_probable_mass(float mag, float m0, float nsat, float *mag_all, float *mass_all, int *indx, int ngal);
float fiber_corrected_galaxy_property(int indx, float *mag_r, float *mag_g, float *prop, int ngal, int *flag);

//local functions
float mass_from_luminosity(float x);
float angular_separation(float a1, float d1, float a2, float d2);
float distance_redshift(float z);
void sdss_luminosity_function(int ngal, float *mag_r, float *redshift, int *indx);
float find_satellites(int i, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, 
		      float x1, int *group_member, int *indx, int ngal, float radius, float mass, 
		      int igrp, float *luminosity, float *nsat_cur, int i1, float *prob_total);
float radial_probability(float mass, float dr, float rad, float ang_rad);
int central_galaxy(int i, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, 
		 float x1, int *group_member, int *indx, int ngal, float radius, float mass, 
		 int igrp, float *luminosity, float *nsat_cur, int i1, float *prob_total);


// global
float GALAXY_DENSITY,
  MSTAR_LIMIT,
  MAGNITUDE,
  REDSHIFT;

int main(int argc, char **argv)
{
  FILE *fp;
  char a,a1[10],a2[10],a3[10],aa[1000];
  int i, j, k,i1, j1, nmatched, idum, n1, n2, ngal, ngal_sample;
  int **galid1, **galid2, xdum;
  float **galpos1, **galpos2, *redshift2, *mass2, *emass2; 
  float *magr, *magg, *ba90, *sersic_n;
  float *d4k, *d4k_err, *con;
  int **ipos1, **ipos2, flag, nrepeat=0;
  float r, rmin, volume, vol_corr;
  float *mag_r, *mag_g, dx, dy, dz, x1, theta, theta_max, rmax, *ra, *dec, 
    *redshift, *mass, *rad, *sigma, *ang_rad, *vmax, *nsat_indi, *ptemp, 
    *luminosity, *group_luminosity;
  int *icollided, *indx, *group_member,*group_index, *reverse_indx, *ka;

  float prob_rad, prob_ang, *prob_total, p, x2, x3, err, mass_tolerance = 0.001, prev_mass, xngrp, maxlum;
  int imass, niter_max=10, igrp, ngrp, *group_center, *temp_group, ngrp_tmp, niter, nsat_tot, imax; 
  double ndens_gal = 0;
  double nsat[160], nhalo[160], nsat_cur, nsub[160]; 
  double nsat2[160], nhalo2[160], mbar2[160], mbars[160]; 
  double mbar[160];

  float *H_delta, *Dn4k, *SFR, *mstellar, *vdisp, *s2n, *sersic, *rexp;
  int BELL_MASSES = 0;

  REDSHIFT = 19200/SPEED_OF_LIGHT;
  MAGNITUDE = -19;
  MSTAR_LIMIT = pow(10.0,9.8);
  if(argc>2)
    {
      MAGNITUDE = atof(argv[1]);
      REDSHIFT = atof(argv[2])/SPEED_OF_LIGHT;
      MSTAR_LIMIT = pow(10.0,atof(argv[3]));
    }
  fprintf(stderr,"MAGNITUDE: %f\nREDSHIFT: %f\nMSTAR: %e\n",MAGNITUDE,REDSHIFT,MSTAR_LIMIT);
  if(argc>4)
    {
      BELL_MASSES = atoi(argv[4]);
      if(BELL_MASSES)fprintf(stderr,"Using Bell et al masses\n");
    }

  volume = 4./3.*PI*(pow(distance_redshift(REDSHIFT),3.0)-pow(distance_redshift(CZMIN/SPEED_OF_LIGHT),3.0))*0.177;
  fprintf(stderr,"Vol %e\n",volume);
  density2halo(0.011);

  for(i=0;i<160;++i)
    nsat[i] = nhalo[i] = mbar[i] = nsub[i] = mbars[i] = 0;

  //read in all the galaxies from the VAGC
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/lss.dr72bright34.dat");
  //fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/bright1/lss.dr72bright1.dat");
  ngal = filesize(fp);
  ra = vector(1,ngal);
  dec = vector(1,ngal);
  redshift = vector(1,ngal);
  sigma = vector(1,ngal);
  ang_rad = vector(1,ngal);
  indx = ivector(1,ngal);
  reverse_indx = ivector(1,ngal);
  prob_total = vector(1,ngal);
  mass = vector(1,ngal);
  rad = vector(1,ngal);
  nsat_indi = vector(1,ngal);
  ptemp = vector(1,ngal);
  luminosity = vector(1,ngal);
  group_center = ivector(1,ngal);
  temp_group = ivector(1,ngal);
  group_member = ivector(1,ngal);
  group_index = ivector(1,ngal);
  group_luminosity = vector(1,ngal);

  for(i=1;i<=ngal;++i) 
    {
      fscanf(fp,"%d %d %d %f %f %f",&idum,&idum,&idum,&ra[i],&dec[i],&redshift[i]);
      ra[i] *= PI/180.;
      dec[i] *= PI/180.;
      indx[i] = i;
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in lss.dat\n");

  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/photo_evolve.dr72bright34.dat");
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
      luminosity[i] = pow(10.0,-0.4*(mag_r[i]+20.44));
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading in photoinfo.dat\n");

  // read in Vmax
  vmax = vector(1,ngal);
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/vmax_evolve.dr72bright34.dat");
  for(i=1;i<=ngal;++i)
    fscanf(fp,"%d %f %f %f",&j,&vmax[i],&x1,&x1);
  fclose(fp);
  fprintf(stderr,"Done reading [%d] lines from vmax file\n",ngal);

  //read in stellar mass
  // it looks like if i just replace luminosity with mstellar, everything works out
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/mstellar.dr72bright34.dat");
  mstellar = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: mstellar file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%d %f",&j,&mstellar[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);

  if(BELL_MASSES)
    {
      for(i=1;i<=ngal;++i)
	mstellar[i] = pow(10.0,-.306 + 1.097*(mag_g[i]-mag_r[i])-0.1 - 0.4*(mag_r[i]-4.64));
    }

  //read in H_deta
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/Hdelta_A.dr72bright34.dat");
  H_delta = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: H_delta file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f",&H_delta[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);

  //read in KA parameters
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/ka.dr72bright34.dat");
  ka = ivector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: KA file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f %f",&x1,&x2);
      ka[i] = 0;
      if(x1>0.2 && x2<pow(10.0,x1)-1)ka[i] = 1;
      fgets(aa,1000,fp);
    }
  fclose(fp);

  //read R_exp
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/galrad.dr72bright34.dat");
  rexp = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: KA file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f %f",&x1,&x2);
      rexp[i] = x1;
      fgets(aa,1000,fp);
    }
  fclose(fp);

  //read in Dn4000
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/dn4k_matched.dat");
  Dn4k = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: Dn4k file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f",&Dn4k[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);

  //read in SSFR
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/sfr.dr72bright34.dat");
  SFR = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: Dn4k file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f",&SFR[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);

  //read in velocity dispersion & signa-to-noise
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/vdisp.dr72bright34.dat");
  vdisp = vector(1,ngal);
  s2n = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: velocity dispersion file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f %f %f",&vdisp[i],&x1,&s2n[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading [%d] lines from vdisp file\n",ngal);

  //read in sersic index
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/sersic_n.dr72bright34.dat");
  sersic = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: mstellar file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%d %f",&j,&sersic[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading [%d] lines from sersic file\n",ngal);

  //read in collision 
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/collided.dr72bright34.dat");
  icollided = ivector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: mstellar file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%d",&icollided[i]);
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading [%d] lines from collided file\n",ngal);

  //read in petrorads 
  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/petrorad.dr72bright34.dat");
  con = vector(1,ngal);
  if(filesize(fp)!=ngal)
    {
      printf("ERROR: mstellar file not correct size: %d %d\n",ngal,filesize(fp));
      exit(0);
    }
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f %f",&x1,&x2);
      con[i] = x2/x1;
      fgets(aa,1000,fp);
    }
  fclose(fp);
  fprintf(stderr,"Done reading [%d] lines from petrorad file\n",ngal);


  fprintf(stderr,"here\n");

  //exit(0);
  for(i=1;i<=ngal;++i)
    {
      if(!(mag_r[i]>MAGNITUDE || redshift[i]>REDSHIFT*SPEED_OF_LIGHT || redshift[i]<CZMIN 
	   || mstellar[i]<MSTAR_LIMIT))
	{
	  x1 = Dn4k[i];
	  x3 = SFR[i];
	  x2 = H_delta[i];
	  if(fabs(Dn4k[i]-1.25)<1.0E-4)
	    x1 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Dn4k,ngal,icollided);
	  if(fabs(Dn4k[i]-1.80)<1.0E-4)
	    x1 = fiber_corrected_galaxy_property(i,mag_r,mag_g,Dn4k,ngal,icollided);
	  if(SFR[i]<-98)
	    x3 = fiber_corrected_galaxy_property(i,mag_r,mag_g,SFR,ngal,icollided);
	  if(fabs(H_delta[i])<1.0E-5)
	    x2 = fiber_corrected_galaxy_property(i,mag_r,mag_g,H_delta,ngal,icollided);

	  printf("GALDAT_CORR %d %f %f %f %f %f %f %e %f %f %f %f %f %f %d %f\n",
		 i,mag_r[i],mag_g[i],redshift[i],x1,x2,x3,
		 mstellar[i],ra[i],dec[i], vdisp[i],s2n[i], sersic[i],con[i],ka[i],rexp[i]);
	  fflush(stdout);
	}
    }
  fflush(stdout);
  fprintf(stderr,"Done with galdata %d\n", ngal);
  exit(0);


  //sort everything by galaxy mass (from most massive to least massive)
  for(i=1;i<=ngal;++i)
    mstellar[i] = -mstellar[i];
  sort2(ngal,mstellar,indx);
  for(i=1;i<=ngal;++i)
    {
      mstellar[i] = -mstellar[i];
      luminosity[i] = mstellar[i];
    }
  fprintf(stderr,"Done sorting\n");

  // put the galaxy magnitudes in the same order
  for(i=1;i<=ngal;++i)
    mag_g[i] = mag_r[i];
  for(i=1;i<=ngal;++i)
    mag_r[i] =  mag_g[indx[i]];

  
  // set up the reverse index list
  fprintf(stderr,"setting up reverse index...\n");
  for(i=1;i<=ngal;++i)
    {
      for(j=1;j<=ngal;++j)
	if(indx[j]==i)break;
      reverse_indx[i] = j;
      printf("%d %d\n",i,j);
    }
  fprintf(stderr,"done setting up reverse index\n");

  // set up the mass/radius of each galaxy
  j = 0;
  k=0;
  if(argc<200)
    {
      for(i1=1;i1<=ngal;++i1)
	{
	  i = indx[i1];
	  ndens_gal += 1/vmax[i];
	  mass[i] = density2halo(ndens_gal);
	  rad[i] = pow(3*mass[i]/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
	  ang_rad[i] = rad[i]/distance_redshift(redshift[i]/SPEED_OF_LIGHT);
	  sigma[i] = sqrt(BIG_G*mass[i]/2.0/rad[i]);
	}
      fprintf(stderr,"Done SHAMming the galaxies\n");
      fprintf(stderr,"Number density: %e\n",ndens_gal);
      fprintf(stderr,"%d in bin 120\n",k);
    }
  else
    {
      // read in the mass from the flux-lim sham file
      fp = openfile(argv[1]);
      for(i=1;i<=ngal;++i)
	{
	  //for(j=1;j<=ngal;++j)
	  // if(indx[j]==i)break;
	  //i1 = j;
	  fscanf(fp,"%s %f %f %f %f",aa,&x1,&x2,&x2,&x3);
	  
	  //if(mag_r[i1]>MAGNITUDE || redshift[indx[i1]]>REDSHIFT*SPEED_OF_LIGHT || redshift[indx[i1]]<CZMIN)continue;
	  //printf("%e %e %e %e %e\n",mass[i],x1,mass[i]/x1, vmax[i],x3);
	  //	  continue;

	  mass[i] = x1;
	  rad[i] = pow(3*mass[i]/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
	  ang_rad[i] = rad[i]/distance_redshift(redshift[i]/SPEED_OF_LIGHT);
	  sigma[i] = sqrt(BIG_G*mass[i]/2.0/rad[i]*(1+redshift[i]/SPEED_OF_LIGHT));
	  //printf("%d %d %e\n",i1,i,mass[i]);
	}
      fclose(fp);
    }

  fprintf(stderr,"here\n");

  j=0;
  for(i1=1;i1<=ngal;++i1)
    {
      group_member[i] = 0;
      i = indx[i1];
      if(mag_r[i1]>MAGNITUDE || redshift[indx[i1]]>(REDSHIFT)*SPEED_OF_LIGHT-CZBUF 
	 || redshift[indx[i1]]<CZMIN+CZBUF || mstellar[i1]<MSTAR_LIMIT)
	prob_total[i] = -1;
      else
	{ prob_total[i] = 0; j++; }
      if(!(mag_r[i1]>MAGNITUDE || redshift[indx[i1]]>REDSHIFT*SPEED_OF_LIGHT 
	   || redshift[indx[i1]]<CZMIN || mstellar[i1]<MSTAR_LIMIT));
    }
  ngal_sample = j;
  fprintf(stderr,"%d gals in sample\n",ngal_sample);
  
  j= 0;
  for(i=1;i<=ngal;++i)
    if(prob_total[i]<0)j++;
  fprintf(stderr,"CHECK: %d w/ negative probability: %d\n",j, ngal-j);

  k=0;
  igrp = 0;
  nsat_tot = 0;
  //go through and find associated galaxies
  for(i1=1;i1<=ngal;++i1)
    { 
      i = indx[i1];
      if(mag_r[i1]>MAGNITUDE)continue;
      if(mstellar[i1]<MSTAR_LIMIT)continue;

      if(i1%10000==0)fprintf(stderr,"%d\n",i1);
      if(mag_r[i1]>MAGNITUDE || redshift[i]/SPEED_OF_LIGHT>REDSHIFT 
	 || redshift[i]<CZMIN || mstellar[i1]<MSTAR_LIMIT)continue;
      
      //if(prob_total[i]>0)continue;
      if(group_member[i])continue;
      //if(prob_total[i]>=1)continue; // see if fully a satellite

      igrp++;
      group_luminosity[igrp] = mstellar[i1];
      group_center[igrp] = i;
      group_member[i] = igrp;

      group_luminosity[igrp] += find_satellites(i,ra,dec,redshift,mag_r,ang_rad[i],sigma[i],group_member,indx,ngal,rad[i],mass[i], igrp, luminosity,&nsat_indi[i],i1,prob_total);

      if(nsat_indi[i]<1)k++;
      nsat_tot+= nsat_indi[i]*(1-prob_total[i]);
      printf("NSAT0 %e %e %f %e\n",mass[i],group_luminosity[igrp],nsat_indi[i],prob_total[i]);
      fflush(stdout);
      
    }
  ngrp = igrp;
  fprintf(stderr,"%d groups, %d n=1 groups, fsat= %.2f \n",ngrp,k,nsat_tot*1./ngal_sample);

  // sort by group luminosity
  for(i=1;i<=ngrp;++i)
    group_luminosity[i] *= -1;
  sort2(ngrp, group_luminosity, group_center);
  for(i=1;i<=ngrp;++i)
    group_luminosity[i] *= -1;

  // now abundance match the rank-ordered total luminosity
  j=0;
  xngrp = 0;
  for(i1=1;i1<=ngrp;++i1)
    {
      i= group_center[i1];

      x1 = mass[i];
      mass[i] = density2host_halo(i1/volume);
      rad[i] = pow(3*mass[i]/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
      ang_rad[i] = rad[i]/distance_redshift(redshift[i]/SPEED_OF_LIGHT);
      sigma[i] = sqrt(BIG_G*mass[i]/2.0/rad[i]*(1+redshift[i]/SPEED_OF_LIGHT));
      
      printf("BOO0 %e %e %f\n",x1,mass[i],nsat_indi[i]);
    }

  // now iterate the groups until convergence
  for(niter=1;niter<=niter_max;++niter)
    {
      //reset the group membership
      j= 0;
      for(i=1;i<=ngal;++i)
	if(prob_total[i]<0)j++;
      //fprintf(stderr,"CHECK: %d w/ negative probability: %d\n",j, ngal-j);
      for(i=1;i<=ngal;++i)
	if(prob_total[i]>=0)prob_total[i]=0;
      for(i=1;i<=ngal;++i)
	group_member[i] = nsat_indi[i] = 0;

      //make a temp list of group centers
      for(i=1;i<=ngrp;++i)
	temp_group[i] = group_center[i];
      ngrp_tmp = ngrp;
      ngrp = 0;
      k=0;
      nsat_tot = 0;

      // find members of the current group list
      for(i1=1;i1<=ngrp_tmp;++i1)
	{
	  i = temp_group[i1];
	  if(group_member[i]>0)continue; //check to see if this gal was claimed as group member
	  ngrp++;
	  group_luminosity[ngrp] = luminosity[reverse_indx[i]];
	  group_member[i] = ngrp;
	  group_center[ngrp] = i;

	  nsat_indi[i] = 0;
	  group_luminosity[ngrp] += find_satellites(i,ra,dec,redshift,mag_r,ang_rad[i],sigma[i],group_member,indx,ngal,rad[i],mass[i], ngrp, luminosity,&nsat_indi[i],i1,prob_total);
	  if(nsat_indi[i]==0)k++;
	  nsat_tot+= nsat_indi[i];

	  printf("NSAT%d %d %d %e %e %f %e\n",niter,ngrp,i,mass[i],group_luminosity[ngrp],nsat_indi[i],prob_total[i]);
	  fflush(stdout);
	}
      
      // look for groups around newly "exposed" galaxies
      
      for(i=1;i<=ngal;++i)
	{
	  if(prob_total[i]<0)continue;
	  if(group_member[i])continue;
	  ngrp++;
	  group_luminosity[ngrp] = luminosity[i];
	  group_member[i] = ngrp;
	  group_center[ngrp] = i;

	  nsat_indi[i] = 0;
	  group_luminosity[ngrp] += find_satellites(i,ra,dec,redshift,mag_r,ang_rad[i],sigma[i],group_member,indx,ngal,rad[i],mass[i], ngrp, luminosity,&nsat_indi[i],i1,prob_total);
	  if(nsat_indi[i]==0)k++;
	  nsat_tot += nsat_indi[i];
	}
      
      fprintf(stderr,"%d groups, %d n=1 groups, fsat= %.2f \n",ngrp,k,nsat_tot*1./ngal_sample);

      for(i=1;i<=ngrp;++i)
	group_index[i] = i;

      // sort by group luminosity
      for(i=1;i<=ngrp;++i)
	group_luminosity[i] *= -1;
      sort2(ngrp, group_luminosity, group_index);
      for(i=1;i<=ngrp;++i)
	group_luminosity[i] *= -1;
      
      // now abundance match the rank-ordered total luminosity
      j=0;
      for(i1=1;i1<=ngrp;++i1)
	{
	  igrp = group_index[i1];
	  i= group_center[igrp];
	  
	  //find new group center
	  if(nsat_indi[i]<-2)
	    {
	      j = central_galaxy(i,ra,dec,redshift,mag_r,ang_rad[i],sigma[i],group_member,indx,ngal,ang_rad[i],mass[i], igrp, luminosity,&nsat_indi[i],i1,prob_total);
	      if(j!=i)
		{
		  group_center[igrp] = j;
		  i = j;
		}
	      printf("CEN%d %d %d %d %f %f\n",niter,igrp,i,j,luminosity[i],luminosity[j]);
	    }

	  x1 = mass[i];
	  mass[i] = density2host_halo(i1/volume);
	  rad[i] = pow(3*mass[i]/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
	  ang_rad[i] = rad[i]/distance_redshift(redshift[i]/SPEED_OF_LIGHT);
	  sigma[i] = sqrt(BIG_G*mass[i]/2.0/rad[i]*(1+redshift[i]/SPEED_OF_LIGHT));
	  
	  //find brightest galaxy in group
	  maxlum = 0;
	  for(j=1;j<=ngal;++j)
	    {
	      if(group_member[j]==igrp)
		if(luminosity[j]>maxlum) 
		  {
		    maxlum = luminosity[j];
		    imax = j;
		  }
	    }

	  printf("GROUP%d %d %d %e %e %f %e %e %e %e %e %e %d %f\n",niter,igrp,i,x1,mass[i],nsat_indi[i],
		 group_luminosity[i1],luminosity[i],prob_total[i],ra[i],dec[i],redshift[i],imax,maxlum);
	  fflush(stdout);

	}

      //print things out in the original order
      //  for(i=1;i<=ngrp;++i)
      //itemp[i] = group_index[i];


      for(i=1;i<=ngal;++i)
	{
	  if(prob_total[i]<0)continue;
	  
	  //for(j=1;j<=ngal;++j)
	  //  if(indx[j]==i)break;
	  //j1 = j;
	  
	  j1 = reverse_indx[i];
	  for(j=1;j<=ngrp;++j)
	    if(group_index[j]==group_member[i])break;
	  igrp = group_index[j];
	  j = group_center[igrp];
	  theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);
	  printf("PROB%d %d %d %d %e %e %e %e %e %e %e %e %e\n",niter,i,igrp,group_center[igrp],mag_r[j1],
		 prob_total[i],mass[i],nsat_indi[i],mass[j],nsat_indi[j],theta/ang_rad[j],theta,ang_rad[j]);
	}


    }
      
  exit(0);

  j=0;
  k=0;
  for(i=1;i<=ngal;++i)
    {
      if(prob_total[i]<0){ k++; continue; }
      if(prob_total[i]>1)j++;
      if(prob_total[i]>1)prob_total[i]=1;

      printf("GAL %e %e %e %e %e\n",mass[i],prob_total[i],nsat_indi[i],vmax[i],nsat_indi[i]);

      imass = (int)(log10(mass[i])*10.0);
      mbar[imass] += mass[i]*(1-prob_total[i]);
      mbars[imass] += mass[i]*prob_total[i];

      nhalo[imass] += (1-prob_total[i]);
      nsub[imass] += prob_total[i];
    }
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


/***************************************
 * estimate the luminosity function using 1/Vmax weighting
 */
void sdss_luminosity_function(int ngal, float *mag_r, float *redshift, int *indx)
{
  int i1,i,j,k,ibin;
  double sum_gal[240], sum_weight[240], sum_vol=0;
  float volume;
  FILE *fp;
  float *vmax,*zmin,*zmax;

  for(i=0;i<240;++i)
    sum_gal[i]=0;

  //test of volume
  fprintf(stderr,"%e %f\n",4./3.*PI*pow(distance_redshift(.1139538)/1000.0,3.0),10.);
  // test is correct
  
  //get the vmax values for each galaxy
  vmax = vector(1,ngal);
  zmin = vector(1,ngal);
  zmax = vector(1,ngal);

  fp = openfile("/Users/tinker/cosmo/SDSS_DATA/DR7_VAGC/vmax_evolve.dr72bright34.dat");
  for(i=1;i<=ngal;++i)
    fscanf(fp,"%d %f %f %f",&j,&vmax[i],&zmin[i],&zmax[i]);
  fclose(fp);
  fprintf(stderr,"Done reading [%d] lines from vmax file\n",ngal);

  // this is the rank-ordered magnitude list  
  for(i1=1;i1<=ngal;++i1)
    {
      i=indx[i1];
      ibin = (int)(fabs(mag_r[i1])*10.0);
      sum_gal[ibin] += 1/vmax[i];
    }
  for(i=0;i<240;++i)
    if(sum_gal[i]>0)printf("%f %e\n",-(i+0.5)/10.0,sum_gal[i]);

  exit(0);
}


/******************************************
 * go through the list of galaxies find satellites for the
 * current central
 */

float find_satellites(int i, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, 
		      float x1, int *group_member, int *indx, int ngal, float radius, float mass, 
		      int igrp, float *luminosity, float *nsat_cur, int i1, float *prob_total)
{
  int j, j1;
  float dx, dy, dz, theta, prob_ang, vol_corr, prob_rad, grp_lum, p0, p0_central;

  *nsat_cur = 0;
  grp_lum = 0;
  for(j1=1;j1<=ngal;++j1)
    {
      j = indx[j1];
      if(j==i)continue;
      if((mag_r[j1]>MAGNITUDE || redshift[j]>REDSHIFT*SPEED_OF_LIGHT || redshift[j]<CZMIN 
	  || luminosity[j1]<MSTAR_LIMIT))continue;
      //if(prob_total[j]<0)
      //printf("NEG %f %f %e\n",mag_r[j1],redshift[j],log10(luminosity[j1]));
      if(prob_total[j]<0)continue;

      //if already in group, skip
      if(group_member[j])continue;

      dx = fabs(ra[i]-ra[j]);
      if(dx>4*theta_max)continue;
      dy = fabs(dec[i]-dec[j]);
      if(dy>4*theta_max)continue;
      dz = fabs(redshift[i] - redshift[j]);
      if(dz>6*x1)continue;
            
      theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);
      
      if(theta>theta_max)continue;
      // FOR MARLA: I want to keep the PDF going out to 2*Rvir
      //if(theta>3*theta_max)continue;
      
      prob_ang = radial_probability(mass,theta,radius,theta_max);
      prob_rad = exp(-dz*dz/(2*x1*x1))*SPEED_OF_LIGHT/(RT2PI*x1);
      
      //fprintf(stdout,"ACK here %d %e %e %e %e %e %e\n",j, dz/x1,theta/theta_max,prob_ang, prob_rad, prob_ang*prob_rad, x1);
      //if(prob_ang*prob_rad<10)continue;

      p0 = (1 - 1/(1+prob_ang*prob_rad/10));
      if(p0<0){
	printf("ZERO %e\n",p0); p0 = 0;
      }
      // marla: if p0 larger than current p0, use this one.
      if(p0>prob_total[j])prob_total[j] = p0;
      // comment out next line
      //prob_total[j] += p0;
      if(prob_total[j]>1)prob_total[j] = 1;

      if(prob_ang*prob_rad<10)continue;
      group_member[j] =igrp;

      grp_lum += luminosity[j1];
      (*nsat_cur)+=1;

      if(igrp==2)printf("BUH %d %f %f %f\n",j,mag_r[j1],redshift[j],CZMIN);

    }
  
  dz = fabs(redshift[i]-CZMIN);
  vol_corr = 1-0.5*erfc(dz/(ROOT2*x1));
  *nsat_cur/= vol_corr;
  grp_lum/=vol_corr;
  dz = fabs(redshift[i]-REDSHIFT);
  vol_corr = 1-0.5*erfc(dz/(ROOT2*x1));
  *nsat_cur/= vol_corr;
  grp_lum/= vol_corr;
  
  return grp_lum;
}

float radial_probability(float mass, float dr, float rad, float ang_rad)
{
  float c, x, rs, delta, f;

  dr = dr*rad/ang_rad;

  c = 10.0*pow(mass/1.0E+14,-0.11);
  rs = rad/c;
  x = dr/rs;

  if(x<1)
    f = 1/(x*x-1)*(1-log((1+sqrt(1-x*x))/x)/(sqrt(1-x*x)));
  if(x==1)
    f = 1.0/3.0;
  if(x>1)
    f = 1/(x*x-1)*(1-atan(sqrt(x*x-1))/sqrt(x*x-1));

  delta = DELTA_HALO/3.0*c*c*c/(log(1+c)-c/(1+c));

  return 1.0/c_on_H0*2*rs*delta*f;
}


int central_galaxy(int i, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, 
		 float x1, int *group_member, int *indx, int ngal, float radius, float mass, 
		 int igrp, float *luminosity, float *nsat_cur, int i1, float *prob_total)
{
  int ii[10000], icnt, j, imax, j1, icen;
  float pot[10000], theta, maxpot, x2;
  
  icen = i;

  // find the galaxies in this group.
  icnt = 0;
  for(i=1;i<=ngal;++i)
    {
      if(group_member[i]==igrp)
	{
	  icnt++;
	  ii[icnt] = i;
	}
    }

  // now do double-sum to get "potential"
  for(i1=1;i1<=icnt;++i1)
    pot[i1] = 0;

  for(i1=1;i1<=icnt;++i1)
    {
      for(j1=i1+1;j1<=icnt;++j1)
	{
	  i = ii[i1];
	  j = ii[j1];
	  x2 = theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);
	  theta = sqrt(theta*theta + radius*radius/400.0);
	  pot[i1] += luminosity[i]*luminosity[j]/theta;
	  pot[j1] += luminosity[i]*luminosity[j]/theta;
	  //if(i1==25 || j1==25)
	  //printf("POT %d %d %e %e %e %e\n",i1,j1,theta,radius/20.0,x2,luminosity[i]*luminosity[j]/theta);
	}
    }

  // find the member with largest potential
  maxpot = 0;
  for(i1=1;i1<=icnt;++i1)
    {
      i = ii[i1];
      //printf("%d %d %d %f %f %f %f %d\n",i1,i,icen,ra[i],dec[i],pot[i1],maxpot,imax);
      if(pot[i1]>maxpot) { maxpot = pot[i1]; imax = i; }
    }

  return imax;

}
