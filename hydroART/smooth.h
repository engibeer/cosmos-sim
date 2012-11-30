/* kd.h */

#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define DARK	1
#define GAS	2
#define STAR	4

typedef struct Particle {
	float r[3];
	float v[3];
	float fMass;
	int iOrder;
	int iMark;
	float fDensity;
	float vMean[3];
	float fVel2;
	} PARTICLE;

typedef struct bndBound {
	float fMin[3];
	float fMax[3];
	} BND;

typedef struct kdNode {
	float fSplit;
	BND bnd;
	int iDim;
	int pLower;
	int pUpper;
	} KDN;

typedef struct kdContext {
	int nBucket;
	int nParticles;
	int nDark;
	int nGas;
	int nStar;
	int bDark;
	int bGas;
	int bStar;
	int nActive;
	float fTime;
	int nLevels;
	int nNodes;
	int nSplit;
	PARTICLE *p;
	KDN *kdNodes;
	int uSecond;
	int uMicro;
	} * KD;


#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	float INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1,INTRSCT_fDist2;\
	INTRSCT_dx = c[cp].bnd.fMin[0]-x;\
	INTRSCT_dx1 = x-c[cp].bnd.fMax[0];\
	INTRSCT_dy = c[cp].bnd.fMin[1]-y;\
	INTRSCT_dy1 = y-c[cp].bnd.fMax[1];\
	INTRSCT_dz = c[cp].bnd.fMin[2]-z;\
	INTRSCT_dz1 = z-c[cp].bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) {\
		INTRSCT_dx1 += lx;\
		if (INTRSCT_dx1 < INTRSCT_dx) {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x+lx;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dx1 > 0.0) {\
		INTRSCT_dx += lx;\
		if (INTRSCT_dx < INTRSCT_dx1) {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x-lx;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		INTRSCT_fDist2 = 0.0;\
		sx = x;\
		}\
	if (INTRSCT_dy > 0.0) {\
		INTRSCT_dy1 += ly;\
		if (INTRSCT_dy1 < INTRSCT_dy) {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y+ly;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dy1 > 0.0) {\
		INTRSCT_dy += ly;\
		if (INTRSCT_dy < INTRSCT_dy1) {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y-ly;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		}\
	if (INTRSCT_dz > 0.0) {\
		INTRSCT_dz1 += lz;\
		if (INTRSCT_dz1 < INTRSCT_dz) {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z+lz;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dz1 > 0.0) {\
		INTRSCT_dz += lz;\
		if (INTRSCT_dz < INTRSCT_dz1) {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z-lz;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		}\
	}


void kdTime(KD,int *,int *);
int kdInit(KD *,int);
int kdReadTipsy(KD,FILE *,int,int,int);
void kdInMark(KD,char *);
void kdBuildTree(KD);
void kdOrder(KD);
void kdFinish(KD);

#endif


/* smooth.h */

#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED
#define RESMOOTH_SAFE  10


typedef struct pqNode {
	float fKey;
	struct pqNode *pqLoser;
	struct pqNode *pqFromInt;
	struct pqNode *pqFromExt;
	struct pqNode *pqWinner;	/* Only used when building initial tree */
	int p;
	float ax;
	float ay;
	float az;
	} PQ;


typedef struct smContext {
	KD kd;
	int nSmooth;
	float fPeriod[3];
	PQ *pq;
	PQ *pqHead;
	float *pfBall2;
	char *iMark;
	int nListSize;
	float *fList;
	int *pList;
	} * SMX;


#define PQ_INIT(pq,n)\
{\
	int PQ_j;\
	for (PQ_j=0;PQ_j<(n);++PQ_j) {\
		if (PQ_j < 2) (pq)[PQ_j].pqFromInt = NULL;\
		else (pq)[PQ_j].pqFromInt = &(pq)[PQ_j>>1];\
		(pq)[PQ_j].pqFromExt = &(pq)[(PQ_j+(n))>>1];\
		}\
	}


#define PQ_BUILD(pq,n,q)\
{\
	int PQ_i,PQ_j;PQ *PQ_t,*PQ_lt;\
	for (PQ_j=(n)-1;PQ_j>0;--PQ_j) {\
		PQ_i = (PQ_j<<1);\
		if (PQ_i < (n)) PQ_t = (pq)[PQ_i].pqWinner;\
		else PQ_t = &(pq)[PQ_i-(n)];\
		++PQ_i;\
		if (PQ_i < (n)) PQ_lt = (pq)[PQ_i].pqWinner;\
		else PQ_lt = &(pq)[PQ_i-(n)];\
		if (PQ_t->fKey < PQ_lt->fKey) {\
			(pq)[PQ_j].pqLoser = PQ_t;\
			(pq)[PQ_j].pqWinner = PQ_lt;\
			}\
		else {\
			(pq)[PQ_j].pqLoser = PQ_lt;\
			(pq)[PQ_j].pqWinner = PQ_t;\
			}\
		}\
	(q) = (pq)[1].pqWinner;\
	}


#define PQ_REPLACE(q)\
{\
	PQ *PQ_t,*PQ_lt;\
	PQ_t = (q)->pqFromExt;\
	while (PQ_t) {\
		if (PQ_t->pqLoser->fKey > (q)->fKey) {\
			PQ_lt = PQ_t->pqLoser;\
			PQ_t->pqLoser = (q);\
			(q) = PQ_lt;\
			}\
		PQ_t = PQ_t->pqFromInt;\
		}\
	}



int smInit(SMX *,KD,int,float *);
void smFinish(SMX);
void smBallSearch(SMX,float,float *);
int  smBallGather(SMX,float,float *);
void smSmooth(SMX,void (*)(SMX,int,int,int *,float *));
void smReSmooth(SMX,void (*)(SMX,int,int,int *,float *));
void smDensity(SMX,int,int,int *,float *);
void smDensitySym(SMX,int,int,int *,float *);
void smMeanVel(SMX,int,int,int *,float *);
void smMeanVelSym(SMX,int,int,int *,float *);
void smVelDisp(SMX,int,int,int *,float *);
void smVelDispSym(SMX,int,int,int *,float *);
void smNull(SMX,int,int,int *,float *);
void smOutDensity(SMX,FILE *);
void smOutMeanVel(SMX,FILE *);
void smOutVelDisp(SMX,FILE *);
void smOutPhase(SMX,FILE *);
void smOutMach(SMX,FILE *);
void smOutSpeed(SMX,FILE *);


/* Daisuke added the routines below */
void Assign_SmoothDen(SMX);
int kdReadHART(KD,int,int,int);

/* link to common blocks in Fortran code */
#define npmax  7000000
extern struct{ double x[npmax]; } part01_;
extern struct{ double y[npmax]; } part02_;
extern struct{ double z[npmax]; } part03_;
extern struct{ double vx[npmax]; } part04_;
extern struct{ double vy[npmax]; } part05_;
extern struct{ double vz[npmax]; } part06_;
extern struct{ float pw[npmax]; } part14_;
extern struct{ int np1; } hpart01_;
extern struct{ 
  float boxh, Om0, Oml0, Omb0, hubble, aexpn, ainit, gamma;
} runparam_;
extern struct {
  float box0, box, pmmsun, fb, rMpc2g, rkpc2g, rg2Mpc, rg2kpc, vg2kms, rg2pMpc, rg2pkpc;
} factors_;

/* dummy density variable in Assign_SmoothDen() */
#define np1max  2000000
extern struct{ float dnb2[np1max]; } hpart04_;

#endif


/* tipsydefs.h */

#define MAXDIM 3
#define forever for(;;)

typedef float Real;

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} ;

struct gas_particle *gas_particles;

struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} ;

struct dark_particle *dark_particles;

struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} ;

struct star_particle *star_particles;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

struct dump header ;


