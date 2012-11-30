#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define KEV_2_ERGS 1.60219e-9	/* convert from keV to ergs */
#define KEV_2_K 1.160485e+7     /* convert from keV to K */
#define DBL_MIN 1.0e-300   /*  convert from keV to K */

double  *Tlam,*Lambda;
int    Tlength;

/* Integrate the RS spectra with the detector response */
double rs_det( double *ex, int *nrs, float *emin, float *de )
{
        int	i,icgs,dlength,j,loc,numt,spec_bins,k,locl,loch,abs_n;
	double	*spec,lambda_cgs,lambda_det,lambda_s,integral,val,cnts2cgs;
	float	*elo,*ehi,*area,spec_emin,spec_de,spec_scale,bandmin,
	  bandmax,emid,te,*spec_e,abs_e[260],abs_sigma[260],y2abs[260],
	  *absorption,sigma,NH;
        char  dirname[300],infile1[300],infile2[300],command[300];
	char  trash[200],junk[200];
	FILE *inp,*pip;
	void  spline(),splint();

	/* H column density [10^20 cm^-2] */
	NH = 0.0;
	icgs = 0;  /* [icgs = 1 (CGS: erg s^-1 cm^-2], 0 [cnts s^-1] */

	*spec = *ex;
	spec_bins = *nrs;
	spec_emin = *emin;
	spec_de = *de;
	
        strcpy(dirname,"/dworkin/home/daisuke/Raymond/EffArea/");
	strcat(infile1,dirname);
	strcat(infile2,dirname);
	strcat(infile1,"chandra_acis-s-bi_0.3-7.0keV.eff");
	strcat(infile2,"cross_sections.dat");
	
	/* calculate length of detector area file */
	sprintf(command,"wc %s",infile1);
	pip=popen(command,"r");
	fscanf(pip,"%d",&dlength);
	dlength--;
	pclose(pip);
	elo=(float *)calloc(dlength,sizeof(float));
	ehi=(float *)calloc(dlength,sizeof(float));
	area=(float *)calloc(dlength,sizeof(float));
	/* read in detector effective area */
	inp=fopen(infile1,"r");
	fgets(trash,200,inp);
	sscanf(trash,"%s %f %f",&junk,&bandmin,&bandmax);
	printf("# %s:  Band is %5.2f:%5.2f with %d energy channels\n",
	  infile1,bandmin,bandmax,dlength);
	for (i=0;i<dlength;i++) 
	  fscanf(inp,"%f %f %f",elo+i,ehi+i,area+i);
	fclose(inp);

	/* read in the galactic absorption cross sections */
	inp=fopen(infile2,"r");
	abs_n=0;while(fgets(trash,200,inp)!=NULL) {
	  if (strncmp(trash,"#",1)) {
	    sscanf(trash,"%f %f",abs_e+abs_n,abs_sigma+abs_n);
	    abs_n++;
	  }
	}
	fclose(inp);
	/* spline this for future use */
	spline(abs_e-1,abs_sigma-1,abs_n,0.0,0.0,y2abs-1);

	/* cycle through rs-spectra, calculating lambdas */
	
	/* exit if spectra don't extend beyond range of detector sensitivity */
	if (bandmin<spec_emin || bandmax>spec_emin+spec_de*spec_bins) {
	  printf("  Detector energy range larger than R-S spectra range\n");
	  exit(0);
	}

	spec=(double *)calloc(spec_bins,sizeof(double));
	spec_e=(float *)calloc(spec_bins,sizeof(float));
	absorption=(float *)calloc(spec_bins,sizeof(float));

	for (j=0;j<spec_bins;j++) {
	  spec_e[j]=spec_emin+(float)j*spec_de;
	  /* calculate the absorption cross section at this energy */
	  if (spec_e[j]<abs_e[0]) sigma=abs_sigma[0];
	  else { 
	    if (spec_e[j]>abs_e[abs_n-1]) sigma=abs_sigma[abs_n-1];
	    else splint(abs_e-1,abs_sigma-1,y2abs-1,abs_n,spec_e[j],&sigma);
	  }
	  absorption[j]=exp(-1.0*sigma*NH);
	  /* RS spectrum is modified by the galactic absorption */ 
	  spec[j]*=absorption[j];
	}

	/* calculate lambda in detector units */
	lambda_det=0.0;
	for (j=0;j<dlength;j++) {
	  /* integrate over spectral emissivity from elo[] to ehi[] */
	  locl=(int)((elo[j]-spec_emin)/spec_de);
	  loch=(int)((ehi[j]-spec_emin)/spec_de);
	  /* detector bin lies within single R-S bin */
	  if (locl==loch) integral=area[j]*(ehi[j]-elo[j])*spec[locl];
	  else {/* detector bin extends over more than one R-S bin */
	    /* take care of the first partial bin */
	    integral=(spec_de-(elo[j]-spec_e[locl]))*spec[locl];
	    /* add all the enclosed bins */
	    for (k=locl+1;k<loch;k++) integral+=spec_de*spec[k];
	    /* take care of the last partial bin */
	    integral+=(ehi[j]-spec_e[loch])*spec[loch];
	    /* now scale by the area of this detector bin */
	    integral*=area[j];
	  }
	  /* convert from cgs to counts */
	  emid=0.5*(elo[j]+ehi[j]);
	  val=integral/(emid*KEV_2_ERGS);
	  lambda_det+=val;
	}
	/* now calculate theoretical lambda in cgs units */
	loc=(bandmin-spec_emin)/spec_de;
	/* initialize using value from first partial bin */
	lambda_cgs=(spec_de-(bandmin-spec_e[loc]))*spec[loc];
	/* add contribution from all the full bins */
	for (j=loc+1;j<spec_bins;j++) {
	  if (bandmax>spec_e[j+1])
	    lambda_cgs+=spec_de*spec[j];
	  else break;
	}
	/* add contribution from last partial bin */
	lambda_cgs+=(bandmax-spec_e[j])*spec[j];
	/* calculate conversion from observed count rate to cgs flux */
	if (lambda_det>0.0) cnts2cgs=lambda_cgs/lambda_det;
	else cnts2cgs=0.0;

	if (icgs ==1) return(lambda_cgs); 
	if (icgs ==0) return(lambda_det);
}


/* Read in Raymond-Smith cooling function for Zmet [Zsun] */
void read_rs()
{
        int i,flag;
	float junk;
        char  infile[300],command[300];
	FILE *inp,*pip;
	
	flag = 1;
        strcpy(infile,"/dworkin/home/daisuke/ARTsim/Data/RS/rs_0.1-2.4_keV_Z0.3.cgs");

        /* Calculate length of cooling functon file */
        sprintf(command,"wc %s",infile); 
        pip=popen(command,"r");
        fscanf(pip,"%d",&Tlength);
        pclose(pip);
	puts(" ");
        printf("# of lines in the input file: %d\n",Tlength);
        /* Allocate momory for Tlamb & Lambda */
        Tlam=(double *)calloc(Tlength,sizeof(double));
        Lambda=(double *)calloc(Tlength,sizeof(double));

        /* Read in Cooling function [erg s^-1 cm^3] */
        printf("Read in Cooling function :  %s \n",&infile); 
        inp=fopen(infile,"r");
        if (flag==0 || flag==1) {
          for (i=0;i<Tlength;i++) fscanf(inp,"%lf %lf \n",&Tlam[i],&Lambda[i]);
	}
        else if (flag==2) {
          for (i=0;i<Tlength;i++) {
            fscanf(inp,"%lf %lf %f %f\n",&Tlam[i],&Lambda[i],&junk,&junk);
	    printf("%lf %lf \n",&Tlam[i],&Lambda[i]);
            Tlam[i]=Tlam[i]*KEV_2_K;
          }
        } 
        fclose(inp);
        puts("Finished reading Cooling function");
        puts(" ");
}

/* Linearly interpolate the Raymond-Smith cooling function   */
/* output value in units of [ne nH * 10**-23 erg/s/cm^3] */ 
void lambda_rs( float *tp, double *lambda )
{ 
        double  val;
        float   t,tmin=1.0e+05;
        int     n;

	t=*tp;

        /* Interpolate in log temperature space */
        n=(int)Tlength*log10(t/Tlam[0])/log10(Tlam[Tlength-1]/Tlam[0]);
        if (n>Tlength-1) {
          printf("Temperature out of range : n=%d t=%7.4f\n",n,t); 
          exit(1); 
        }
        else if (n<0 || t<tmin)  
          val=DBL_MIN;  /* Cold gas does not emit X-ray much : DBL_MIN=2.225074e-
308 */
        else if (n==0) 
          val=Lambda[0]+(t-Tlam[0])*(Lambda[1]-Lambda[0])/(Tlam[1]-Tlam[0]);
        else if (n==Tlength-1) 
          val=Lambda[Tlength-1];
        else { /* interpolate in log space */
          val=Lambda[n]+(t-Tlam[n])*(Lambda[n+1]-Lambda[n])/(Tlam[n+1]-Tlam[n]);
          /* printf("Check %e %le \n",t,val); */ 
        }
	val = val * pow(10.0,23);
	*lambda=val;
}


/* Read in Oxygen emission lines (OVII & OVIII) */
void load_oxygen(int *iflag)
{
        int i, ioxygen;
	static double temp[30] = { 5.0, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 
				   5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 
				   6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 
				   7.5, 7.6, 7.7, 7.8, 7.9, 8.0 };

	static double OVII_561[30] = { 3.5408e-31, 1.4620e-30, 8.9348e-31,
				       4.7585e-31, 3.3952e-31, 1.1450e-29,
				       3.3500e-28, 4.5681e-27, 3.5454e-26,
				       1.7582e-25, 5.9681e-25, 1.3709e-24,
				       2.0075e-24, 1.7497e-24, 9.4474e-25,
				       3.9948e-25, 1.6123e-25, 6.7380e-26,
				       2.9675e-26, 1.3767e-26, 6.6797e-27,
				       3.3714e-27, 1.7598e-27, 9.4206e-28,
				       5.1726e-28, 2.8995e-28, 1.6552e-28,
				       9.6301e-29, 5.7081e-29, 3.4543e-29 };

	static double OVII_568[30] = { 1.1060e-31, 4.4174e-31, 2.6521e-31, 
				       1.6450e-31, 5.5556e-30, 1.2341e-28,
				       1.1166e-27, 6.3205e-27, 2.6655e-26,
				       8.9324e-26, 2.3625e-25, 4.5907e-25,
				       5.9646e-25, 4.7562e-25, 2.4032e-25,
				       9.6718e-26, 3.7636e-26, 1.5313e-26,
				       6.6103e-27, 3.0195e-27, 1.4443e-27,
				       7.2054e-28, 3.7168e-28, 1.9638e-28,
				       1.0622e-28, 5.8498e-29, 3.2699e-29,
				       1.8559e-29, 1.0685e-29, 6.2519e-30 };
	
	static double OVII_574[30] = { 9.3913e-32, 3.6485e-31, 2.1605e-31,
				       1.1213e-31, 6.4130e-31, 3.4950e-29,
				       7.3558e-28, 8.2586e-27, 5.7568e-26,
 				       2.7167e-25, 9.0484e-25, 2.0637e-24,
				       3.0059e-24, 2.6024e-24, 1.3924e-24,
				       5.8148e-25, 2.3101e-25, 9.4759e-26,
				       4.0864e-26, 1.8543e-26, 8.7997e-27, 
				       4.3398e-27, 2.2141e-27, 1.1587e-27,
				       6.2203e-28, 3.4084e-28, 1.9007e-28,
				       1.0792e-28, 6.2334e-29, 3.6688e-29 };

	static double OVII_665[30] = { 1.4746e-32, 6.1955e-32, 3.8149e-32, 
				       2.0455e-32, 1.5586e-32, 8.2178e-31,
				       2.9022e-29, 4.6804e-28, 4.2337e-27,
				       2.4245e-26, 9.3517e-26, 2.3821e-25,
				       3.7669e-25, 3.4680e-25, 1.9428e-25,
				       8.4018e-26, 3.4291e-26, 1.4367e-26,
				       6.3049e-27, 2.9038e-27, 1.3968e-27,
				       6.9775e-28, 3.6048e-28, 1.9111e-28,
				       1.0396e-28, 5.7771e-29, 3.2700e-29,
				       1.8862e-29, 1.1077e-29, 6.6342e-30 };

	static double OVIII_653[30] = { 3.8250e-31, 1.4900e-30, 8.8533e-31,
				       	4.5866e-31, 2.0969e-31, 8.0000e-33,
				       	9.8752e-35, 7.3396e-32, 1.4535e-29,
				       	9.9705e-28, 2.7051e-26, 3.0294e-25,
				       	1.4794e-24, 3.3308e-24, 3.9066e-24,
				       	3.2012e-24, 2.3137e-24, 1.6326e-24,
				       	1.1623e-24, 8.4245e-25, 6.2225e-25,
				       	4.6722e-25, 3.5584e-25, 2.7380e-25,
				       	2.1266e-25, 1.6633e-25, 1.3080e-25,
				       	1.0341e-25, 8.2172e-26, 6.5693e-26 };     

	ioxygen = *iflag;

        /* Allocate momory for Tlamb & Lambda */
	Tlength=30;
        Tlam=(double *)calloc(Tlength,sizeof(double));
        Lambda=(double *)calloc(Tlength,sizeof(double));

        /* Read in Cooling function [erg s^-1 cm^3] */
        printf("Loading oxygen emissivities...\n");
	for (i=0;i<Tlength;i++) {
	  Tlam[i]=pow(10.0,temp[i]);
	  if ( ioxygen == 1 ) Lambda[i]=OVII_561[i];
	  if ( ioxygen == 2 ) Lambda[i]=OVII_568[i];
	  if ( ioxygen == 3 ) Lambda[i]=OVII_574[i];
	  if ( ioxygen == 4 ) Lambda[i]=OVII_665[i];
	  if ( ioxygen == 5 ) Lambda[i]=OVIII_653[i];
	  /* printf("%le %le \n",Tlam[i],Lambda[i]); */
	}
        /* puts("Finished loadig the oxygen emissivities"); */
        puts(" ");
}


/* Linearly interpolate the Oxygen line emissions */
/* output value in units of [ne**2 erg s^-1 cm^-3] */ 
void lambda_oxygen( float *tp, double *lambda )
{ 
        double  val;
        float   t,tmin=1.0e+05;
        int     n;

	t=*tp;

        /* Interpolate in log temperature space */
        n=(int)Tlength*log10(t/Tlam[0])/log10(Tlam[Tlength-1]/Tlam[0])-1;
        if (n<0 || t<tmin || n>Tlength-1)  
          val=0.0;
        else if (n==0) 
          val=Lambda[0]+(t-Tlam[0])*(Lambda[1]-Lambda[0])/(Tlam[1]-Tlam[0]);
        else if (n==Tlength-1) 
          val=Lambda[Tlength-1];
        else { /* interpolate in log space */
          val=Lambda[n]+(t-Tlam[n])*(Lambda[n+1]-Lambda[n])/(Tlam[n+1]-Tlam[n]);
          /* printf("Check %e %le \n",t,val); */ 
        }
	
	/* printf("%d %d %e %e %e\n",Tlength,n,t,Tlam[0],Tlam[Tlength-1]);
	   printf("%d %e %e %e %e %e %e\n",n,t,Tlam[n],Tlam[n+1],Lambda[n],Lambda[n+1],val); */
	
	*lambda=val;	
}

void spline(x,y,n,yp1,ypn,y2)
float x[],y[],yp1,ypn,y2[];
int n;
{
	int i,k;
	float p,qn,sig,un,*u,*vector();
	void free_vector();

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}

void splint(xa,ya,y2a,n,x,y)
float xa[],ya[],y2a[],x,*y;
int n;
{
	int klo,khi,k;
	float h,b,a;
	void nrerror();

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

/* for different function name convention */
double rs_det_(double *ex, int *nrs, float *emin, float *de) { rs_det(ex,nrs,emin,de); }
void read_rs_() { read_rs(); }
void lambda_rs_(float *tp, double *lambda) { lambda_rs(tp,lambda); }
void load_oxygen_(int *iflag) { load_oxygen(iflag); }
void lambda_oxygen_(float *tp, double *lambda) { lambda_oxygen(tp,lambda); }


