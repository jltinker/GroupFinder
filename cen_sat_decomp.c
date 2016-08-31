#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define OMEGA_M 0.3
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

#define CZMIN 7000
#define CZBUF 0
#define REDSHIFT (19200.0/SPEED_OF_LIGHT)
#define MAGNITUDE -19.0

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
float halo_abundance(float m);
float subhalo_abundance(float m);
float subhalo_abundance2(float m);
float func_nhalo(float m);
float func_subhalo1(float mhost);
void scatter_test(void);
float most_probable_mass(float mag, float m0, float nsat, float *mag_all, float *mass_all, int *indx, int ngal);

//local functions
float mass_from_luminosity(float x);
float angular_separation(float a1, float d1, float a2, float d2);
float distance_redshift(float z);
void sdss_luminosity_function(int ngal, float *mag_r, float *redshift, int *indx);
float find_satellites(int i, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, 
		      float x1, float *prob_total, int *indx, int ngal);

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
  float r, rmin, volume, vol_corr;
  float *mag_r, *mag_g, dx, dy, dz, x1, theta, theta_max, rmax, *ra, *dec, *redshift, *mass, *rad, *sigma, *ang_rad, *vmax, *nsat_indi, *ptemp;
  int *icollided, *indx;

  float prob_rad, prob_ang, *prob_total, p, x2, x3, err, mass_tolerance = 0.001, prev_mass;
  int imass, niter_max=100; 
  double ndens_gal = 0;
  double nsat[160], nhalo[160], nsat_cur, nsub[160]; 
  double nsat2[160], nhalo2[160], mbar2[160], mbars[160]; 
  double mbar[160];

  //scatter_test();
  //exit(0);

  //REDSHIFT = 17000/SPEED_OF_LIGHT;

  // flux-limited sample
  //REDSHIFT = 1;
  //CZMIN = 0;
  // MAGNITUDE = -4;

  volume = 4./3.*PI*(pow(distance_redshift(REDSHIFT),3.0)-pow(distance_redshift(CZMIN/SPEED_OF_LIGHT),3.0))*0.19312;
  fprintf(stderr,"Vol %e\n",volume);
  volume = 4./3.*PI*(pow(distance_redshift(31000./SPEED_OF_LIGHT),3.0)-pow(distance_redshift(CZMIN/SPEED_OF_LIGHT),3.0))*0.19312;
  fprintf(stderr,"Vol %e\n",volume);
  volume = 4./3.*PI*(pow(distance_redshift(48000./SPEED_OF_LIGHT),3.0)-pow(distance_redshift(CZMIN/SPEED_OF_LIGHT),3.0))*0.19312;
  fprintf(stderr,"Vol %e\n",volume);
  volume = 4./3.*PI*(pow(distance_redshift(73000./SPEED_OF_LIGHT),3.0)-pow(distance_redshift(CZMIN/SPEED_OF_LIGHT),3.0))*0.19312;
  fprintf(stderr,"Vol %e\n",volume);
  density2halo(0.011);

  /*
  for(i=100;i<161;++i)
    {
      fprintf(stdout,"%e %e %e %e\n",pow(10.0,(i+0.5)/10),
	      halo_abundance(pow(10.0,(i+0.5)/10.0))*pow(10.0,(i+0.5)/10),
	      subhalo_abundance2(log(pow(10.0,(i+0.5)/10.0))),
	      qromo(subhalo_abundance2,log(pow(10.0,i/10.0)),(i+1.0)/10.0*log(10),midpnt)/0.1/log(10));
      fflush(stdout);
    }
  exit(0);
  for(i=100;i<161;++i)
    {
      fprintf(stdout,"%e %e %e %e\n",pow(10.0,(i+0.5)/10),
	      qromo(halo_abundance,pow(10.0,i/10.0),pow(10.0,(i+1.0)/10.0),midpnt)*volume,
		    qromo(func_nhalo,log(pow(10.0,i/10.0)),(i+1.0)/10.0*log(10),midpnt)*volume,
		    qromo(subhalo_abundance2,log(pow(10.0,i/10.0)),(i+1.0)/10.0*log(10),midpnt)*volume);
      fflush(stdout);
    }
  //exit(0);
  */

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
  prob_total = vector(1,ngal);
  mass = vector(1,ngal);
  rad = vector(1,ngal);
  nsat_indi = vector(1,ngal);
  ptemp = vector(1,ngal);

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
      //mag_r[i] = mag_r[i] + Q0*(1+Q1*(redshift[i]/SPEED_OF_LIGHT-QZ0))*(redshift[i]/SPEED_OF_LIGHT-QZ0);
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


  //sort everything by luminosity
  sort2(ngal,mag_r,indx);
  fprintf(stderr,"Done sorting\n");

  // see what you get for the luminosity function
  //sdss_luminosity_function(ngal,  mag_r,  redshift, indx);

  // set up the mass/radius of each galaxy
  j = 0;
  k=0;
  if(argc<2)
    {
      for(i1=1;i1<=ngal;++i1)
	{
	  i = indx[i1];
	  if(mag_r[i1]>MAGNITUDE || redshift[indx[i1]]>REDSHIFT*SPEED_OF_LIGHT || redshift[indx[i1]]<CZMIN)continue;
	  j++;
	  //mass[i] = mass_from_luminosity(mag_r[i]);
	  //mass[i] = density2halo(j/volume);
	  //if(log10(mass[i])>=12.0 && log10(mass[i])<12.1) { k++; fprintf(stderr,"%e %d\n",mass[i],(int)(log10(mass[i])*10.0)); }
	  ndens_gal += 1/vmax[i];
	  mass[i] = density2halo(ndens_gal);
	  rad[i] = pow(3*mass[i]/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
	  ang_rad[i] = rad[i]/distance_redshift(redshift[i]/SPEED_OF_LIGHT);
	  sigma[i] = sqrt(BIG_G*mass[i]/2.0/rad[i]);
	  //printf("BOO %f %e %e %e %e %e %d %e\n",mag_r[i1],mass[i],rad[i],ang_rad[i],sigma[i],distance_redshift(redshift[i]/SPEED_OF_LIGHT),j,func_nhalo(log(mass[i]))*0.1*log(10)*volume);
	  //qromo(func_nhalo,log(mass[i]),log(1.0E+16),midpnt)*volume);
	  
	}
      fprintf(stderr,"Done SHAMming the galaxies\n");
      fprintf(stderr,"Number density: %d %e\n",j,j/volume);
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


  j=0;
  for(i1=1;i1<=ngal;++i1)
    {
      i = indx[i1];
      if(mag_r[i1]>MAGNITUDE || redshift[indx[i1]]>(REDSHIFT)*SPEED_OF_LIGHT-CZBUF || redshift[indx[i1]]<CZMIN+CZBUF)
	prob_total[i] = -1;
      else
	prob_total[i] = 0;

      //if(ra[i]<140*ANG || ra[i]>230*ANG || dec[i]<0 || dec[i]>45*ANG)prob_total[i] = -1;

      if(mag_r[i1]<=MAGNITUDE && redshift[indx[i1]]<(REDSHIFT)*SPEED_OF_LIGHT-CZBUF && redshift[indx[i1]]>CZMIN+CZBUF)j++;
	//if(!(ra[i]<140*ANG || ra[i]>230*ANG || dec[i]<0 || dec[i]>45*ANG))j++;
    }
  fprintf(stderr,"%d gals in sample\n",j);

  k=0;
  //go through and find associated galaxies
  for(i1=1;i1<=ngal;++i1)
    {
      i = indx[i1];

      i = 486355;
      for(j=1;j<=ngal;++j)
	if(indx[j]==i)break;
      i1 = j;

      if(mag_r[i1]>MAGNITUDE)break;

      if(i1%10000==0)fprintf(stderr,"%d\n",i1);
      if(mag_r[i1]>MAGNITUDE || redshift[i]/SPEED_OF_LIGHT>REDSHIFT || redshift[i]<CZMIN)continue;

      imass = (int)(log10(mass[i])*10.0);

      x1 = sqrt(sigma[i]*sigma[i] + 100.*100.*rad[i]*rad[i]);
      nsat_cur = 0;

      theta_max = ang_rad[i];//*pow(2.0*3.0,THIRD);

      //check if cur gal has too high sat prob
      if(prob_total[i]>1)prob_total[i]=1;

      // now iterate to convergence of the halo mass
      prev_mass = mass[i];
      for(j=1;j<=niter_max;++j)
	{
	  x1 = sqrt(sigma[i]*sigma[i] + 100.*100.*rad[i]*rad[i]);
	  printf("PRE%d %d %e %e %e %e\n",i1,j,mass[i],rad[i],ang_rad[i],x1);
	  nsat_cur = find_satellites(i1,ra,dec,redshift,mag_r,ang_rad[i],x1,ptemp,indx,ngal);
	  printf("POST %e\n",nsat_cur);
	  if(nsat_cur==0)break;
	  mass[i] = most_probable_mass(mag_r[i1],mass[i],nsat_cur,mag_r,mass,indx,ngal);

	  if(isnan(mass[i]))
	    {
	      //printf("BOO NAN\n");
	      mass[i] = prev_mass;
	      break;
	    }

	  rad[i] = pow(3*mass[i]/(4.*PI*DELTA_HALO*RHO_CRIT*OMEGA_M),1.0/3.0);
	  ang_rad[i] = rad[i]/distance_redshift(redshift[i]/SPEED_OF_LIGHT);
	  sigma[i] = sqrt(BIG_G*mass[i]/2.0/rad[i]);

	  err = fabs(mass[i] - prev_mass)/prev_mass;
	  printf("BOO%d %d %e %e %e %e\n",i1,j,mass[i],prev_mass,nsat_cur,err);
	  fflush(stdout);
	  prev_mass = mass[i];

	  if(err<mass_tolerance)break;
	}

      //now assign the probabilities to the sat halos
      for(j1=i1+1;j1<=ngal;++j1)
	{
	  j = indx[j1];
	  prob_total[j] += ptemp[j];
	}
      
      nsat_indi[i] = nsat_cur;
      imass = (int)(log10(mass[i])*10.0);
      nsat[imass]+=nsat_cur;

      printf("BOO %e %e %e\n",mass[i],nsat_cur,prob_total[i]);
      exit(0);
    }

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

      //nhalo[imass] += (1-prob_total[i]);/vmax[i];
      //nsub[imass] += prob_total[i]/vmax[i];
    }

  /*
  for(i=0;i<160;++i)
    nsat2[i] = nhalo2[i] = 0;
  for(i=0;i<160;++i)
    {      
      for(j=i-12;j<=i+12;++j)
	{
	  if(j<0)continue;
	  if(j>=160)continue;
	  p = G0*exp(-(j*1.0-i*1.0)*(j-i*1.0)/10.0/10.0/(2*0.3*0.3))*0.3*2.302;
	  //printf("%d %d %e\n",j,i,p);
	  //mbar2[j] += mbar[i]*(j-1)/p;
	  nsat2[j] += nsat[i]*p;
	  nhalo2[j] += nhalo[i]*p;
	}
      mbar2[i] = pow(10.0,(i+0.5)/10.0);
    }
  */
  // exit(0);

  fprintf(stderr,"%d %d\n",j,ngal-k);
  for(i=0;i<160;++i)
    if(nhalo[i]>0)
      //printf("%e %e %e %e %e %d\n",mbar2[i],1.0,nsat2[i]/nhalo2[i],nhalo2[i],nsat2[i],i);
      printf("%e %e %e %e %e %e %d\n",mbar[i]/nhalo[i],1.0,nsat[i]/nhalo[i],nhalo[i],nsat[i],nsub[i],i);

  exit(0);

  fprintf(stderr,"%d %d\n",j,ngal-k);
  for(i=0;i<160;++i)
    if(nhalo[i]>0)
      //printf("%e %e %e %e %e %d\n",mbar2[i],1.0,nsat2[i]/nhalo2[i],nhalo2[i],nsat2[i],i);
      printf("%e %e %e %e %e %d %e %e %e %e\n",mbar[i]/nhalo[i],1.0,nsat[i]/nhalo[i],nhalo[i],nsat[i],i,
	     halo_abundance(mbar[i]/nhalo[i])*mbar[i]/nhalo[i]*0.1*log(10)*volume,
	     nsub[i],subhalo_abundance((mbar[i]/nhalo[i]))*mbar[i]/nhalo[i]*0.1*log(10)*volume,
	     qromo(subhalo_abundance2,log(pow(10.0,i/10.0)),(i+1.0)/10.0*log(10),midpnt)*volume);
  //	     qromo(subhalo_abundance2,log(pow(10.0,i/10.0)),log(1.0E+15),midpnt)*volume);

  //mbars[i]/nsub[i],subhalo_abundance((mbars[i]/nsub[i]))*mbars[i]/nsub[i]*0.1*log(10)*volume

  
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
      fp = openfile("massfunc.wmap1");
      //fp = openfile("sham_lum2mass.dat");
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
  fprintf(stderr,"%e %f\n",4./3.*PI*pow(distance_redshift(.1139538)/1000.0,3.0),10);
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

float find_satellites(int i1, float *ra, float *dec, float *redshift, float *mag_r, float theta_max, 
		      float x1, float *prob_total, int *indx, int ngal)
{
  int i, j, j1;
  float dx, dy, dz, theta, prob_ang, nsat_cur, vol_corr, prob_rad;

  i = indx[i1];

  nsat_cur = 0;
  for(j1=i1+1;j1<=ngal;++j1)
    {
      j = indx[j1];
      if(redshift[j]>REDSHIFT*SPEED_OF_LIGHT-CZBUF || mag_r[j1]>MAGNITUDE || redshift[i]<CZMIN+CZBUF)continue;
      
      dx = fabs(ra[i]-ra[j]);
      if(dx>2*theta_max)continue;
      dy = fabs(dec[i]-dec[j]);
      if(dy>2*theta_max)continue;
      dz = fabs(redshift[i] - redshift[j]);
      if(dz>6*x1)continue;
      
      
      theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);
      if(theta>theta_max)continue;
      
      prob_ang = radial_probability(mass,theta,radius,theta_max);
      //prob_ang *= (1-prob_total[i]);
      prob_rad = exp(-dz*dz/(2*x1*x1));
      nsat_cur += prob_ang*prob_rad;
      prob_total[j] = prob_ang*prob_rad;
      printf("GUH %d %e\n",j,nsat_cur);
    }
  dz = fabs(redshift[i]-CZMIN);
  vol_corr = 1-0.5*erfc(dz/(ROOT2*x1));
  nsat_cur/= vol_corr;
  dz = fabs(redshift[i]-REDSHIFT);
  vol_corr = 1-0.5*erfc(dz/(ROOT2*x1));
  nsat_cur/= vol_corr;

  return nsat_cur;
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

  delta = 200.0/3.0*c*c*c/(log(1+c)-c/(1+c));

  return 1.0/c_on_H0*2*rs*delta*f;
}
