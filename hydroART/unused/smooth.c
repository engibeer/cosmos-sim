/* ***************************************************************** */
/* smooth.c                                                          */
/*    Original code was modified to work with HART code              */
/*       SmoothDen() : a driver to compute smoothed density          */
/*       kdReadHART  : read position and velocities of DM particles  */
/*       Assign_SmoothDen : assign smoothed DM density to dnb2       */
/*                                                                   */
/*    Calling Fortran routines/variables :                           */
/*      On the NCSA, do not add _ after the routines/variable names  */
/*      On the local machine, add _ after the routines/variable names*/
/*                                                                   */
/*    Last modified (Daisuke 10/27/03)                               */ 
/* ***************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include "smooth.h"

#define BIGNOTQUITEMAXFLOAT  ((float)1.0e+37)
#define MAX_ROOT_ITTR	32


/* main.c */
void usage(void)
{
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"smooth [-s <nSmooth>[dgs]] [-b <nBucket>] [-g]\n");
	fprintf(stderr,"   [-o <Output Name>] [-p <xyzPeriod>]\n");
	fprintf(stderr,"   [-px <xPeriod>] [-py <yPeriod>] [-pz <zPeriod>]\n");
	fprintf(stderr,"   [-do <MarkFile>]\n");
	fprintf(stderr,"   [density] [meanvel] [speed] [veldisp] [mach]\n");
	fprintf(stderr,"   [phase] [all] [null]\n\n");
	fprintf(stderr,"Input taken from stdin in tipsy binary format.\n");
	fprintf(stderr,"SEE MAN PAGE: smooth(1) for more information.\n");
	exit(1);
	}

void SmoothDen( int *nsm, int *ngrid )
{
	KD kd;
	SMX smx;
	int nBucket,nSmooth,i,j;
	FILE *fp;
	char ach[256],achFile[80],achMark[80];
	float fPeriod[3];
	int bDensity,bMeanVel,bVelDisp,bPhase,bMach,bSpeed,bNull,bSym;
	int bDark,bGas,bStar;
	int bMark;
	char *p,*q;

	int argcc;
	char argvv[20][200];
	int ic, i1, i2, i3, i4, ismooth;
	extern void write_density_hf( int *ismooth );

	/* Daisuke added lines below */
	argcc = 6;
	strcpy(argvv[1],"-s");
	sprintf(argvv[2],"%d%s",*nsm,"d"); /* use nsm part. to compute local density */
	strcpy(argvv[3],"-p");
	sprintf(argvv[4],"%d",*ngrid);
	strcpy(argvv[5],"density");
	
	printf("smooth %s %s %s %s %s \n",argvv[1],argvv[2],
	       argvv[3],argvv[4],argvv[5]);
	
   	nBucket = 16;
	nSmooth = 64;
	bDensity = 0;
	bMeanVel = 0;
	bVelDisp = 0;
	bSpeed = 0;
	bMach = 0;
	bPhase = 0;
	bNull = 0;
	bSym = 1;
	bDark = 1;
	bGas = 1;
	bStar = 1;
	bMark = 0;
	
	/* Name of the output file */
	ic = (int)(runparam_.aexpn*1000);
	i1 = (int)(ic / 1000); 
	i2 = (int)((ic - (ic/1000)*1000)/100); 
	i3 = (int)((ic - (ic/100)*100) / 10); 
	i4 = (int)((ic - (ic/10)*10));
	ismooth = 1;
	strcpy(achFile,"rho_smooth_a");
	sprintf(achFile,"%s%d%s%d%d%d",achFile,i1,".",i2,i3,i4);
        
	i = 1;
	for (j=0;j<3;++j) fPeriod[j] = BIGNOTQUITEMAXFLOAT;
	while (i < argcc) {
		if (!strcmp(argvv[i],"-b")) {
			++i;
			if (i >= argcc) usage();
			nBucket = atoi(argvv[i]);
			++i;
			}
		else if (!strcmp(argvv[i],"-s")) {
			++i;
			if (i >= argcc) usage();
			p = argvv[i];
			while (isdigit(*p)) ++p;
			q = p;
			if (isalpha(*p)) {
				bDark = 0;
				bGas = 0;
				bStar = 0;
				}
			while (isalpha(*p)) {
				switch (*p) {
				case 'd':
					bDark = 1;
					break;
				case 'g':
					bGas = 1;
					break;
				case 's':
					bStar = 1;
					break;
				default:
					usage();
					}
				++p;
				}
			*q = 0;
			nSmooth = atoi(argvv[i]);
			++i;
			}
		else if (!strcmp(argvv[i],"-o")) {
			++i;
			if (i >= argcc) usage();
			strcpy(achFile,argvv[i]);
			++i;
			}
		else if (!strcmp(argvv[i],"-do")) {
			++i;
			if (i >= argcc) usage();
			strcpy(achMark,argvv[i]);
			bMark = 1;
			bSym = 0;	/* Symmetrical kernal is inconsistent here! */
			++i;
			}
		else if (!strcmp(argvv[i],"-g")) {
			bSym = 0;
			++i;
			}
		else if (!strcmp(argvv[i],"-p")) {
			++i;
			if (i >= argcc) usage();
			fPeriod[0] = atof(argvv[i]);
			fPeriod[1] = atof(argvv[i]);
			fPeriod[2] = atof(argvv[i]);
			++i;
			}
		else if (!strcmp(argvv[i],"-px")) {
			++i;
			if (i >= argcc) usage();
			fPeriod[0] = atof(argvv[i]);
			++i;
			}
		else if (!strcmp(argvv[i],"-py")) {
			++i;
			if (i >= argcc) usage();
			fPeriod[1] = atof(argvv[i]);
			++i;
			}
		else if (!strcmp(argvv[i],"-pz")) {
			++i;
			if (i >= argcc) usage();
		    fPeriod[2] = atof(argvv[i]);
			++i;
			}
		else if (!strcmp(argvv[i],"density")) {
			bDensity |= 3;
		    ++i;
			}
		else if (!strcmp(argvv[i],"meanvel")) {
			bDensity |= 1;
			bMeanVel |= 3;
			++i;
			}
		else if (!strcmp(argvv[i],"veldisp")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bVelDisp |= 3;
			++i;
			}
		else if (!strcmp(argvv[i],"phase")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bVelDisp |= 1;
			bPhase |= 2;
			++i;
			}
		else if (!strcmp(argvv[i],"mach")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bVelDisp |= 1;
			bMach |= 2;
			++i;
			}
		else if (!strcmp(argvv[i],"speed")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bSpeed |= 2;
			++i;
			}
		else if (!strcmp(argvv[i],"null")) {
			bNull |= 1;
			++i;
			}
		else if (!strcmp(argvv[i],"all")) {
			bDensity |= 3;
			bMeanVel |= 3;
			bVelDisp |= 3;
			bPhase |= 2;
			bMach |= 2;
			bSpeed |= 2;
			++i;
			}
		else usage();
		}
	kdInit(&kd,nBucket);
	/* kdReadTipsy(kd,stdin,bDark,bGas,bStar); */
	kdReadHART(kd,bDark,bGas,bStar);
	if (bMark) kdInMark(kd,achMark);
	kdBuildTree(kd);
	smInit(&smx,kd,nSmooth,fPeriod);
	if (bNull&1) {
		smSmooth(smx,smNull);
		smReSmooth(smx,smNull);
		}
	if (bSym) {
		if (bDensity&1) smSmooth(smx,smDensitySym);
		if (bMeanVel&1) smReSmooth(smx,smMeanVelSym);
		if (bVelDisp&1) smReSmooth(smx,smVelDispSym);
		}
	else {
		if (bDensity&1) smSmooth(smx,smDensity);
		if (bMeanVel&1) smReSmooth(smx,smMeanVel);
		if (bVelDisp&1) smReSmooth(smx,smVelDisp);
		}
	kdOrder(kd);
	if (bDensity&2) {
		strcpy(ach,achFile);
		strcat(ach,".den");
		puts(" Assigning density to dnb2");
		Assign_SmoothDen(smx);
		printf(" Writing : %s \n",ach);
		/* fp = fopen(ach,"w"); 
		   assert(fp != NULL); 
		   smOutDensity(smx,fp);
		   fclose(fp); */
		/* use Fortran routine for I/O below */
		/* write_density_hf( &ismooth ); */
		}
	if (bMeanVel&2) {
		strcpy(ach,achFile);
		strcat(ach,".mvl");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutMeanVel(smx,fp);
		fclose(fp);
		}
	if (bSpeed&2) {
		strcpy(ach,achFile);
		strcat(ach,".spd");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutSpeed(smx,fp);
		fclose(fp);
		}
	if (bVelDisp&2) {
		strcpy(ach,achFile);
		strcat(ach,".dsp");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutVelDisp(smx,fp);
		fclose(fp);
		}
	if (bMach&2) {
		strcpy(ach,achFile);
		strcat(ach,".mch");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutMach(smx,fp);
		fclose(fp);
		}
	if (bPhase&2) {
		strcpy(ach,achFile);
		strcat(ach,".phs");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutPhase(smx,fp);
		fclose(fp);
		}
	smFinish(smx);
	kdFinish(kd);
	puts("End of SmoothDen()");
	return;
	}

/* Daisuke modified KdReadTipsy to read in HART variables */
int kdReadHART(KD kd,int bDark,int bGas,int bStar)
{
	int i,j,nCnt;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;

	/* fread(&h,sizeof(struct dump),1,fp); */
	h.ndim  = 3;
	h.time  = (double)runparam_.aexpn;
	h.nbodies = (int)hpart01_.np1;
	h.nsph  = 0;
	h.nstar = 0;
	h.ndark = h.nbodies - h.nstar - h.nsph;
	printf(" Correcting for fb = %f \n",factors_.fb);
	/* printf("%lf %d %d %d %d \n",h.time,h.nbodies,h.nsph,h.nstar,h.ndark); */

	kd->nParticles = h.nbodies;
	kd->nDark = h.ndark;
	kd->nGas = h.nsph;
	kd->nStar = h.nstar;
	kd->fTime = h.time;
	kd->nActive = 0;
	if (bDark) kd->nActive += kd->nDark;
	if (bGas) kd->nActive += kd->nGas;
	if (bStar) kd->nActive += kd->nStar;
	kd->bDark = bDark;
	kd->bGas = bGas;
	kd->bStar = bStar;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<h.nsph;++i) {
	        /* fread(&gp,sizeof(struct gas_particle),1,fp); */
		if (bGas) {
			kd->p[nCnt].fMass = gp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = gp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = gp.vel[j];
			++nCnt;
			}
		}
	for (i=0;i<h.ndark;++i) {
	        /* fread(&dp,sizeof(struct dark_particle),1,fp); */
		if (bDark) {
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			kd->p[nCnt].fMass = part14_.pw[nCnt] / (1.0-factors_.fb); 
		        /* kd->p[nCnt].fMass = part14_.pw[nCnt]; */
			kd->p[nCnt].r[0] = (float)part01_.x[nCnt];
			kd->p[nCnt].r[1] = (float)part02_.y[nCnt];
			kd->p[nCnt].r[2] = (float)part03_.z[nCnt];
			kd->p[nCnt].v[0] = (float)part04_.vx[nCnt];
			kd->p[nCnt].v[1] = (float)part05_.vy[nCnt];
			kd->p[nCnt].v[2] = (float)part06_.vz[nCnt];
			/* if ( nCnt < 5 ) {
			  printf("%d %f %f %f %f %f %f %f \n",
				 nCnt,kd->p[nCnt].fMass,
				 kd->p[nCnt].r[0],kd->p[nCnt].r[1],
				 kd->p[nCnt].r[2],kd->p[nCnt].v[0],
				 kd->p[nCnt].v[1],kd->p[nCnt].v[2]);
				 } */
			++nCnt;
			}
		}
	for (i=0;i<h.nstar;++i) {
	        /* fread(&sp,sizeof(struct star_particle),1,fp); */
		if (bStar) {
			kd->p[nCnt].fMass = sp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = sp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = sp.vel[j];
			++nCnt;
			}
		}
	return(kd->nParticles);
	}


/* Assign smoothed DM density to dnb2 for halo finding */
void Assign_SmoothDen(SMX smx)
{
	int i,iCnt;

	iCnt = 0;
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
		        if (smx->kd->p[iCnt].iMark) 
			  hpart04_.dnb2[iCnt] = smx->kd->p[iCnt].fDensity;
			else 
			  hpart04_.dnb2[iCnt] = 0.0;
			/* if ( iCnt < 5 ) {
			  printf("as : %d %f \n",iCnt,hpart04_.dnb2[iCnt]);
			  } */
			++iCnt;
			}
		else {
		  printf(" bDark = ",smx->kd->bDark);
		  exit(0);
		  }
	        }
	}

	
/* kd.c */
void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket)
{
	KD kd;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	kd->p = NULL;
	kd->kdNodes = NULL;
	*pkd = kd;
	return(1);
	}


int kdReadTipsy(KD kd,FILE *fp,int bDark,int bGas,int bStar)
{
	int i,j,nCnt;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;

	fread(&h,sizeof(struct dump),1,fp);
	kd->nParticles = h.nbodies;
	kd->nDark = h.ndark;
	kd->nGas = h.nsph;
	kd->nStar = h.nstar;
	kd->fTime = h.time;
	kd->nActive = 0;
	if (bDark) kd->nActive += kd->nDark;
	if (bGas) kd->nActive += kd->nGas;
	if (bStar) kd->nActive += kd->nStar;
	kd->bDark = bDark;
	kd->bGas = bGas;
	kd->bStar = bStar;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<h.nsph;++i) {
		fread(&gp,sizeof(struct gas_particle),1,fp);
		if (bGas) {
			kd->p[nCnt].fMass = gp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = gp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = gp.vel[j];
			++nCnt;
			}
		}
	for (i=0;i<h.ndark;++i) {
		fread(&dp,sizeof(struct dark_particle),1,fp);
		if (bDark) {
			kd->p[nCnt].fMass = dp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = dp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = dp.vel[j];
			++nCnt;
			}
		}
	for (i=0;i<h.nstar;++i) {
		fread(&sp,sizeof(struct star_particle),1,fp);
		if (bStar) {
			kd->p[nCnt].fMass = sp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = sp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = sp.vel[j];
			++nCnt;
			}
		}
	return(kd->nParticles);
	}


void kdInMark(KD kd,char *pszFile)
{
	FILE *fp;
	char ach[80];
	int i,iCnt,iDum;

	fp = fopen(pszFile,"r");
	if (!fp) {
		fprintf(stderr,"Could not open mark array, %s\n",pszFile);
		exit(1);
		}
	fgets(ach,80,fp);	/* ignore the array header! */
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	fclose(fp);
	}


void kdSelect(KD kd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	double v;
	int i,j;

	p = kd->p;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void kdUpPass(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		kdUpPass(kd,l);
		kdUpPass(kd,u);
		kdCombine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->p[u].r[j];
			c[iCell].bnd.fMax[j] = kd->p[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->p[pj].r[j];
				if (kd->p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->p[pj].r[j];
				}
			}
		}
	}


void kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes != NULL) free(kd->kdNodes);
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->p[0].r[j];
		bnd.fMax[j] = kd->p[0].r[j];
		}
	for (i=1;i<kd->nActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->p[i].r[j]) 
				bnd.fMin[j] = kd->p[i].r[j];
			else if (bnd.fMax[j] < kd->p[i].r[j])
				bnd.fMax[j] = kd->p[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelect(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->p[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	kdUpPass(kd,ROOT);
	}


int cmpParticles(const void *v1,const void *v2)
{
	PARTICLE *p1=(PARTICLE *)v1,*p2=(PARTICLE *)v2;
	
	return(p1->iOrder - p2->iOrder);
	}


void kdOrder(KD kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLE),cmpParticles);
	}


void kdFinish(KD kd)
{
	free(kd->p);
	free(kd->kdNodes);
	free(kd);
	}


/* smooth.c */
int smInit(SMX *psmx,KD kd,int nSmooth,float *fPeriod)
{
	SMX smx;
	KDN *root;
	int pi,j;
	int bError=0;

	root = &kd->kdNodes[ROOT];
	assert(root != NULL);
	/*
	 ** Check to make sure that the bounds of the simulation agree 
	 ** with the period specified, if not cause an error.
	 */
	for (j=0;j<3;++j) {
		if (root->bnd.fMax[j] - root->bnd.fMin[j] > fPeriod[j]) {
			fprintf(stderr,"ERROR(smInit):Bounds of the simulation volume exceed\n");
			fprintf(stderr,"exceed the period specified in the %c-dimension.\n",'x'+j);
			bError = 1;
			}
		}
	if (bError) exit(1);
	assert(nSmooth <= kd->nActive);
	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
	smx->kd = kd;
	smx->nSmooth = nSmooth;
	smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
	assert(smx->pq != NULL);
	PQ_INIT(smx->pq,nSmooth);
	smx->pfBall2 = (float *)malloc((kd->nActive+1)*sizeof(int));
	assert(smx->pfBall2 != NULL);
	smx->iMark = (char *)malloc(kd->nActive*sizeof(char));
	assert(smx->iMark);
	smx->nListSize = smx->nSmooth+RESMOOTH_SAFE;
	smx->fList = (float *)malloc(smx->nListSize*sizeof(float));
	assert(smx->fList != NULL);
	smx->pList = (int *)malloc(smx->nListSize*sizeof(int));
	assert(smx->pList != NULL);
	/*
	 ** Set for Periodic Boundary Conditions.
	 */
	for (j=0;j<3;++j) smx->fPeriod[j] = fPeriod[j];
	/*
	 ** Initialize arrays for calculated quantities.
	 */
	for (pi=0;pi<smx->kd->nActive;++pi) {
		smx->kd->p[pi].fDensity = 0.0;
		for (j=0;j<3;++j) smx->kd->p[pi].vMean[j] = 0.0;
		smx->kd->p[pi].fVel2 = 0.0;
		}
	*psmx = smx;	
	return(1);
	}


void smFinish(SMX smx)
{
	free(smx->pfBall2);
	free(smx->iMark);
	free(smx->pq);
	free(smx->fList);
	free(smx->pList);
	free(smx);
	}


void smBallSearch(SMX smx,float fBall2,float *ri)
{
	KDN *c;
	PARTICLE *p;
	int cell,cp,ct,pj;
	float fDist2,dx,dy,dz,lx,ly,lz,sx,sy,sz,x,y,z;
	PQ *pq;

	c = smx->kd->kdNodes;
	p = smx->kd->p;
	pq = smx->pqHead;
	x = ri[0];
	y = ri[1];
	z = ri[2];
	lx = smx->fPeriod[0];
	ly = smx->fPeriod[1];
	lz = smx->fPeriod[2];
	cell = ROOT;
	/*
	 ** First find the "local" Bucket.
	 ** This could mearly be the closest bucket to ri[3].
	 */
	while (cell < smx->kd->nSplit) {
		if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	/*
	 ** Now start the search from the bucket given by cell!
	 */
	for (pj=c[cell].pLower;pj<=c[cell].pUpper;++pj) {
		dx = x - p[pj].r[0];
		dy = y - p[pj].r[1];
		dz = z - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < fBall2) {
			if (smx->iMark[pj]) continue;
			smx->iMark[pq->p] = 0;
			smx->iMark[pj] = 1;
			pq->fKey = fDist2;
			pq->p = pj;
			pq->ax = 0.0;
			pq->ay = 0.0;
			pq->az = 0.0;
			PQ_REPLACE(pq);
			fBall2 = pq->fKey;
			}
		}
	while (cell != ROOT) {
		cp = SIBLING(cell);
		ct = cp;
		SETNEXT(ct);
		while (1) {
			INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz);
			/*
			 ** We have an intersection to test.
			 */
			if (cp < smx->kd->nSplit) {
				cp = LOWER(cp);
				continue;
				}
			else {
				for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
					dx = sx - p[pj].r[0];
					dy = sy - p[pj].r[1];
					dz = sz - p[pj].r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						if (smx->iMark[pj]) continue;
						smx->iMark[pq->p] = 0;
						smx->iMark[pj] = 1;
						pq->fKey = fDist2;
						pq->p = pj;
						pq->ax = sx - x;
						pq->ay = sy - y;
						pq->az = sz - z;
						PQ_REPLACE(pq);
						fBall2 = pq->fKey;
						}
					}
				}
		GetNextCell:
			SETNEXT(cp);
			if (cp == ct) break;
			}
		cell = PARENT(cell);
		}
	smx->pqHead = pq;
	}


int smBallGather(SMX smx,float fBall2,float *ri)
{
	KDN *c;
	PARTICLE *p;
	int pj,nCnt,cp,nSplit;
	float dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;

	c = smx->kd->kdNodes;
	p = smx->kd->p;
	nSplit = smx->kd->nSplit;
	lx = smx->fPeriod[0];
	ly = smx->fPeriod[1];
	lz = smx->fPeriod[2];
	x = ri[0];
	y = ri[1];
	z = ri[2];
	nCnt = 0;
	cp = ROOT;
	while (1) {
		INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz);
		/*
		 ** We have an intersection to test.
		 */
		if (cp < nSplit) {
			cp = LOWER(cp);
			continue;
			}
		else {
			for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
				dx = sx - p[pj].r[0];
				dy = sy - p[pj].r[1];
				dz = sz - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 < fBall2) {
					smx->fList[nCnt] = fDist2;
					smx->pList[nCnt++] = pj;
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}


void smSmooth(SMX smx,void (*fncSmooth)(SMX,int,int,int *,float *))
{
	KDN *c;
	PARTICLE *p;
    PQ *pq,*pqLast;
	int cell;
	int pi,pin,pj,pNext,nCnt,nSmooth;
	float dx,dy,dz,x,y,z,h2,ax,ay,az;


	for (pi=0;pi<smx->kd->nActive;++pi) {
		if (smx->kd->p[pi].iMark) smx->pfBall2[pi] = -1.0;
		else smx->pfBall2[pi] = 1.0;	/* pretend it is already done! */
		}
	smx->pfBall2[smx->kd->nActive] = -1.0; /* stop condition */
	for (pi=0;pi<smx->kd->nActive;++pi) {
		smx->iMark[pi] = 0;
		}
	pqLast = &smx->pq[smx->nSmooth-1];
	c = smx->kd->kdNodes;
	p = smx->kd->p;
	nSmooth = smx->nSmooth;
	/*
	 ** Initialize Priority Queue.
	 */
	pin = 0;
	pNext = 1;
	ax = 0.0;
	ay = 0.0;
	az = 0.0;
	for (pq=smx->pq,pj=0;pq<=pqLast;++pq,++pj) {
		smx->iMark[pj] = 1;
		pq->p = pj;
		pq->ax = ax;
		pq->ay = ay;
		pq->az = az;
		}
	while (1) {
		if (smx->pfBall2[pin] >= 0) {
			/*
			 ** Find next particle which is not done, and load the
			 ** priority queue with nSmooth number of particles.
			 */
			while (smx->pfBall2[pNext] >= 0) ++pNext;
			/*
			 ** Check if we are really finished.
			 */
			if (pNext == smx->kd->nActive) break;
			pi = pNext;
			++pNext;
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			/*
			 ** First find the "local" Bucket.
			 ** This could mearly be the closest bucket to ri[3].
			 */
			cell = ROOT;
			while (cell < smx->kd->nSplit) {
				if (p[pi].r[c[cell].iDim] < c[cell].fSplit)
					cell = LOWER(cell);
				else
					cell = UPPER(cell);
				}
			/*
			 ** Remove everything from the queue.
			 */
			smx->pqHead = NULL;
			for (pq=smx->pq;pq<=pqLast;++pq) smx->iMark[pq->p] = 0;
			/*
			 ** Add everything from pj up to and including pj+nSmooth-1.
			 */
			pj = c[cell].pLower;
			if (pj > smx->kd->nActive - nSmooth)
				pj = smx->kd->nActive - nSmooth;
			for (pq=smx->pq;pq<=pqLast;++pq) {
				smx->iMark[pj] = 1;
				dx = x - p[pj].r[0];
				dy = y - p[pj].r[1];
				dz = z - p[pj].r[2];
				pq->fKey = dx*dx + dy*dy + dz*dz;
				pq->p = pj++;
				pq->ax = 0.0;
				pq->ay = 0.0;
				pq->az = 0.0;
				}
			PQ_BUILD(smx->pq,nSmooth,smx->pqHead);
			}
		else {
			/*
			 ** Calculate the priority queue using the previous particles!
			 */
			pi = pin;
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			smx->pqHead = NULL;
			for (pq=smx->pq;pq<=pqLast;++pq) {
				pq->ax -= ax;
				pq->ay -= ay;
				pq->az -= az;
				dx = x + pq->ax - p[pq->p].r[0];
				dy = y + pq->ay - p[pq->p].r[1];
				dz = z + pq->az - p[pq->p].r[2];
				pq->fKey = dx*dx + dy*dy + dz*dz;
				}
			PQ_BUILD(smx->pq,nSmooth,smx->pqHead);
			ax = 0.0;
			ay = 0.0;
			az = 0.0;
			}
		smBallSearch(smx,smx->pqHead->fKey,p[pi].r);
		smx->pfBall2[pi] = smx->pqHead->fKey;
		/*
		 ** Pick next particle, 'pin'.
		 ** Create fList and pList for function 'fncSmooth'.
		 */
		pin = pi;
		nCnt = 0;
		h2 = smx->pqHead->fKey;
		for (pq=smx->pq;pq<=pqLast;++pq) {
			if (pq == smx->pqHead) continue;
			smx->pList[nCnt] = pq->p;
			smx->fList[nCnt++] = pq->fKey;
			if (smx->pfBall2[pq->p] >= 0) continue;
			if (pq->fKey < h2) {
				pin = pq->p;
				h2 = pq->fKey;
				ax = pq->ax;
				ay = pq->ay;
				az = pq->az;
				}
			}
		(*fncSmooth)(smx,pi,nCnt,smx->pList,smx->fList);
		}
	}


void smReSmooth(SMX smx,void (*fncSmooth)(SMX,int,int,int *,float *))
{
	PARTICLE *p;
	int pi,nSmooth;

	p = smx->kd->p;
	for (pi=0;pi<smx->kd->nActive;++pi) {
		if (p[pi].iMark == 0) continue;
		/*
		 ** Do a Ball Gather at the radius of the most distant particle
		 ** which is smDensity sets in smx->pBall[pi].
		 */
		nSmooth = smBallGather(smx,smx->pfBall2[pi],p[pi].r);
		(*fncSmooth)(smx,pi,nSmooth,smx->pList,smx->fList);
		}
 	}


void smDensity(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float ih2,r2,rs,fDensity;
	int i,pj;

	ih2 = 4.0/smx->pfBall2[pi];
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		pj = pList[i];
		r2 = fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		fDensity += rs*smx->kd->p[pj].fMass;
		}
	smx->kd->p[pi].fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	}


void smDensitySym(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs;
	int i,pj;

	ih2 = 4.0/smx->pfBall2[pi];
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		pj = pList[i];
		r2 = fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		smx->kd->p[pi].fDensity += rs*smx->kd->p[pj].fMass;
		smx->kd->p[pj].fDensity += rs*smx->kd->p[pi].fMass;
		}
	}


void smDensityScat(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs;
	int i,pj;

	ih2 = 4.0/smx->pfBall2[pi];
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		pj = pList[i];
		r2 = fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		smx->kd->p[pj].fDensity += rs*smx->kd->p[pi].fMass;
		}
	}


void smMeanVel(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs;
	int i,j,pj;

	ih2 = 4.0/smx->pfBall2[pi];
	fNorm = M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		pj = smx->pList[i];
		r2 = smx->fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		for (j=0;j<3;++j) {
			smx->kd->p[pi].vMean[j] += rs*smx->kd->p[pj].fMass/
				smx->kd->p[pj].fDensity*smx->kd->p[pj].v[j];
			}
		}
	}


void smMeanVelSym(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs;
	int i,j,pj;

	ih2 = 4.0/smx->pfBall2[pi];
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		pj = smx->pList[i];
		r2 = smx->fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		for (j=0;j<3;++j) {
			smx->kd->p[pi].vMean[j] += rs*smx->kd->p[pj].fMass/
				smx->kd->p[pj].fDensity*smx->kd->p[pj].v[j];
			smx->kd->p[pj].vMean[j] += rs*smx->kd->p[pi].fMass/
				smx->kd->p[pi].fDensity*smx->kd->p[pi].v[j];
			}
		}
	}


void smVelDisp(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs,tv2;
	int i,j,pj;

	ih2 = 4.0/smx->pfBall2[pi];
	fNorm = M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		pj = smx->pList[i];
		r2 = smx->fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		tv2 = 0.0;
		for (j=0;j<3;++j) {
			tv2 += (smx->kd->p[pj].v[j] - smx->kd->p[pi].vMean[j])*
				(smx->kd->p[pj].v[j] - smx->kd->p[pi].vMean[j]);
			}
		smx->kd->p[pi].fVel2 += rs*smx->kd->p[pj].fMass/
			smx->kd->p[pj].fDensity*tv2;
		}
	}	


void smVelDispSym(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	float fNorm,ih2,r2,rs,tv2;
	int i,j,pj;

	ih2 = 4.0/smx->pfBall2[pi];
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		pj = smx->pList[i];
		r2 = smx->fList[i]*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		tv2 = 0.0;
		for (j=0;j<3;++j) {
			tv2 += (smx->kd->p[pj].v[j] - smx->kd->p[pi].vMean[j])*
				(smx->kd->p[pj].v[j] - smx->kd->p[pi].vMean[j]);
			}
		smx->kd->p[pi].fVel2 += rs*smx->kd->p[pj].fMass/
			smx->kd->p[pj].fDensity*tv2;
		tv2 = 0.0;
		for (j=0;j<3;++j) {
			tv2 += (smx->kd->p[pi].v[j] - smx->kd->p[pj].vMean[j])*
				(smx->kd->p[pi].v[j] - smx->kd->p[pj].vMean[j]);
			}
		smx->kd->p[pj].fVel2 += rs*smx->kd->p[pi].fMass/
			smx->kd->p[pi].fDensity*tv2;
		}
	}


void smOutDensity(SMX smx,FILE *fp)
{
	int i,iCnt;

	fprintf(fp,"%d\n",smx->kd->nParticles);
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].fDensity);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].fDensity);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].fDensity);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	}


void smOutMeanVel(SMX smx,FILE *fp)
{
	int i,iCnt;

	fprintf(fp,"%d\n",smx->kd->nParticles);
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[0]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[0]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[0]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[1]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[1]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[1]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[2]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[2]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark)
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].vMean[2]);
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	}


void smOutVelDisp(SMX smx,FILE *fp)
{
	int i,iCnt;

	fprintf(fp,"%d\n",smx->kd->nParticles);
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark) {
				fprintf(fp,"%.8g\n",sqrt(smx->kd->p[iCnt].fVel2));
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark) {
				fprintf(fp,"%.8g\n",sqrt(smx->kd->p[iCnt].fVel2));
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark) {
				fprintf(fp,"%.8g\n",sqrt(smx->kd->p[iCnt].fVel2));
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	}


void smOutPhase(SMX smx,FILE *fp)
{
	int i,iCnt;

	fprintf(fp,"%d\n",smx->kd->nParticles);
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark) {
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].fDensity/
						sqrt(smx->kd->p[iCnt].fVel2));
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark) {
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].fDensity/
						sqrt(smx->kd->p[iCnt].fVel2));
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark) {
				fprintf(fp,"%.8g\n",smx->kd->p[iCnt].fDensity/
						sqrt(smx->kd->p[iCnt].fVel2));
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	}


void smOutMach(SMX smx,FILE *fp)
{
	int i,iCnt;
	float vx,vy,vz,v2,mach;

	fprintf(fp,"%d\n",smx->kd->nParticles);
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark) {
				vx = smx->kd->p[iCnt].vMean[0];
				vy = smx->kd->p[iCnt].vMean[1];
				vz = smx->kd->p[iCnt].vMean[2];
				v2 = vx*vx + vy*vy + vz*vz;
				mach = sqrt(v2/smx->kd->p[iCnt].fVel2);
				fprintf(fp,"%.8g\n",mach);
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark) {
				vx = smx->kd->p[iCnt].vMean[0];
				vy = smx->kd->p[iCnt].vMean[1];
				vz = smx->kd->p[iCnt].vMean[2];
				v2 = vx*vx + vy*vy + vz*vz;
				mach = sqrt(v2/smx->kd->p[iCnt].fVel2);
				fprintf(fp,"%.8g\n",mach);
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark) {
				vx = smx->kd->p[iCnt].vMean[0];
				vy = smx->kd->p[iCnt].vMean[1];
				vz = smx->kd->p[iCnt].vMean[2];
				v2 = vx*vx + vy*vy + vz*vz;
				mach = sqrt(v2/smx->kd->p[iCnt].fVel2);
				fprintf(fp,"%.8g\n",mach);
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	}


void smOutSpeed(SMX smx,FILE *fp)
{
	int i,iCnt;
	float vx,vy,vz,v2,speed;

	fprintf(fp,"%d\n",smx->kd->nParticles);
	iCnt = 0;
	for (i=0;i<smx->kd->nGas;++i) {
		if (smx->kd->bGas) {
			if (smx->kd->p[iCnt].iMark) {
				vx = smx->kd->p[iCnt].vMean[0];
				vy = smx->kd->p[iCnt].vMean[1];
				vz = smx->kd->p[iCnt].vMean[2];
				v2 = vx*vx + vy*vy + vz*vz;
				speed = sqrt(v2);
				fprintf(fp,"%.8g\n",speed);
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nDark;++i) {
		if (smx->kd->bDark) {
			if (smx->kd->p[iCnt].iMark) {
				vx = smx->kd->p[iCnt].vMean[0];
				vy = smx->kd->p[iCnt].vMean[1];
				vz = smx->kd->p[iCnt].vMean[2];
				v2 = vx*vx + vy*vy + vz*vz;
				speed = sqrt(v2);
				fprintf(fp,"%.8g\n",speed);
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	for (i=0;i<smx->kd->nStar;++i) {
		if (smx->kd->bStar) {
			if (smx->kd->p[iCnt].iMark) {
				vx = smx->kd->p[iCnt].vMean[0];
				vy = smx->kd->p[iCnt].vMean[1];
				vz = smx->kd->p[iCnt].vMean[2];
				v2 = vx*vx + vy*vy + vz*vz;
				speed = sqrt(v2);
				fprintf(fp,"%.8g\n",speed);
				}
			else fprintf(fp,"0\n");
			++iCnt;
			}
		else fprintf(fp,"0\n");
		}
	}


void smNull(SMX smx,int pi,int nSmooth,int *pList,float *fList)
{
	return;
	}


/* Fortran-C wrapper */
/* Fortran calls should be in lower case! */
/* for different function name convention */
void smoothden_fortran_( int *nsm, int *ngrid ) {  SmoothDen( nsm, ngrid );  } 
