#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "array.h"
#include "complex.h"
#include "mie.h"
#include "nrutil.h"
#include <time.h>

#define ALIVE      	1
#define DEAD       	0
#define	NN      	100                  /*sg---> number of loops running(pixels 100*100)--------------------->*/ 
#define THRESHOLD   0.01                /* used in roulette */
#define CHANCE      0.1                 /* used in roulette */

#define RandomNum (double) RandomGen(1, 0, NULL)
#define SIGN(x) ((x)>=0 ? 1:-1)
#define InitRandomGen (double) RandomGen(0, 1, NULL)      

/* Declare Subroutines */


/*----comment by me (sandeep) will start with sg-------------------------->*/ 

/* sg---> Rotates the Stokes vector S by an angle phi, storing the result in S2.--------------------->*/ 
void rotSphi(double* S, double phi, double* S2);
double	RandomGen(char Type, long Seed, long *Status);

/* sg---> Multiplies the Stokes vector S by a scalar theta, storing the result in S2.*/ 
void multS(double* S, double theta, double* S2);

/*sg---> Rotates the plane defined by XX, YY, and ZZ by an angle phi and theta, storing the result in XX2, YY2, and ZZ2.*/ 
void rotateXXYY(double* XX, double* YY,double* ZZ, double phi, double theta, double* XX2, double* YY2,double* ZZ2);

/*sg---> Updates the Stokes vector U by an angle phi and theta, storing the result in U2.*/ 
void updateU(double* U, double phi, double theta, double* U2);

/*sg---> Computes the sine and cosine of an angle.*/ 

double sincos(double *x);

/*************** MAIN ****************************************/	
int main() {
	double pi = 3.1415926535897932384;
	
	/* Mie theory stuff */

	/*sg---> Radius of the sphere , wavelength, scattering.--------------------->*/
	double radius,lambda, A;
	/*sg---> Number of angles, i.--------------------->*/
	long nangles,i;
	/*sg---> Stokes vector.--------------------->*/
	struct complex m;
	/*sg---> Complex numbers.--------------------->*/
	struct complex*s1=NULL;
	/*sg---> Stokes vector.--------------------->*/
	struct complex*s2=NULL;



	/*sg---> Mie scattering parameters.--------------------->*/
	double *mu=NULL;
	/*sg---> Scattering parameters.--------------------->*/
	double x,qext,qsca,qback,g, rho, vol,dy, dx, hw;
	double nre_p, nim_p, nre_med, nim_med;
	FILE *target;

	/*sg---> Photon position.--------------------->*/
	double jjj;
	
	/* E field stuff */
	double	phi, theta,I,I0;
	int 	ithedeg;
	double	IT, QT, UT, VT;
	double	IR_1, QR_1, UR_1, VR_1;

	double	**IR, **QR, **UR, **VR;
    
	/* Propagation parameters */
	double	y,z;		/* photon position.  x already declared. Also, incrementals & max range. */
	double	s;          /* step sizes. s = -log(RND)/mus [cm] */
	long	i_photon;   /* current photon */
	long	Nphotons;   /* number of photons in simulation */
	short   photon_status;  /*  = ALIVE=1 or DEAD=0 */

	/* other variables */
	double	mua;        /* absorption coefficient [cm^-1] */
	double	mus;        /* scattering coefficient [cm^-1] */
	double	musp;       /* reduced scattering coefficient [cm^-1] */
	double  albedo;     /* albedo of tissue */

	/* dummy variables */
	double  rnd;        /* assigned random value 0-1 */
	int     MM;
	double  W,absorb;  /* photon weight */
	double  slabsize; 
	int 	j,ix,iy;
	double  cos22,sin22,costheta,sini,cosi;
	
/**** allocate matrices and arrays *******/
	double	*U, *U2;
	double	*S;     	/* */
	double	*S0;     	/* */
	double	*S2;     	/* */
	double	*s11=NULL;
	double	*s12=NULL;
	double	*s33=NULL;
	double	*s43=NULL;
	double	*IQUV;     	/* [I, Q, U, V] Stokes Vector */
	double	 start_time,finish_time,temp;

start_time = clock();

/*sg---> Initialize Stokes vectors and complex numbers--------------------->*/
MM = NN - 1;
/*sg---> Allocate memory for the Stokes vectors--------------------->*/
U      = new_darray(3);
U2     = new_darray(3);
S      = new_darray(4);
S0     = new_darray(4);
S2     = new_darray(4);/* dummy S*/


IQUV   = new_darray(4);                   /*sg---> Stokes vector[I, Q, U, V]--------------------->*/
IR     = dmatrix(0, MM, 0, MM); /* [0:MM] */
QR     = dmatrix(0, MM, 0, MM); /* [0:MM] */
UR     = dmatrix(0, MM, 0, MM); /* [0:MM] */
VR     = dmatrix(0, MM, 0, MM); /* [0:MM] */




/**** end  allocate matrices and arrays *******/


/* CHOOSE MIE SCATTERING parameters */
radius  	= 2.03/2; /* microns */
lambda 		= 0.355; /* microns */
rho 		= 1.152e-4;/*Dilution 1*/
Nphotons	= 1e6;
mua 		= 0.0; /*�a  */

/* ------------------------*/
nre_p   	= 1.59;
nim_p  	 	= 0;
nre_med		= 1.33;
nim_med 	= 0.0;
nangles 	= 1000;                                /*sk-- angle defined--------------------->*/


/* Setup MIE SCATTERING parameters */

mu  = new_darray(nangles);
s1  = new_carray(nangles);
s2  = new_carray(nangles);
s11 = new_darray(nangles);
s12 = new_darray(nangles);
s33 = new_darray(nangles);
s43 = new_darray(nangles);


m.re = nre_p/nre_med;
m.im = 0.0;
x    = 2*pi*radius/(lambda/nre_med);                         /*sk---> from Mie theory--------------------->*/
vol  = 4.0/3*pi*radius*radius*radius;
A    = pi*radius*radius;

for(i=0;i<=nangles;i++)
	mu[i] = cos(pi*i/nangles);
	s11=new_darray(nangles);
	s12=new_darray(nangles);
	s33=new_darray(nangles);
	s43=new_darray(nangles);		
	s1=new_carray(nangles);
	s2=new_carray(nangles);
	
	Mie(x,m,mu,nangles,s1,s2,&qext,&qsca,&qback,&g); /* <---- Call Mie program ----- */

	mus 	= qsca*A*rho*1e4; /* Mus is in cm^-1 */
	musp 	= mus*(1-g);/* [cm^-1] */
	albedo 	= mus/(mus + mua);
	free_darray(mu);

	printf("Polarized Monte Carlo\n dia=%5.5f;\n mus=%5.5f;\n g=%5.5f;\n  rho=%5.5f;\n",radius*2,mus,g,rho); 
	
	slabsize	= 4/mus;
	
/*Scattering parameters s11 s12 s33 s43*/

for(i=0;i<=nangles;++i){
	s11[i] = 0.5*cabbs(s2[i])*cabbs(s2[i]) + 0.5*cabbs(s1[i])*cabbs(s1[i]);
	s12[i] = 0.5*cabbs(s2[i])*cabbs(s2[i]) - 0.5*cabbs(s1[i])*cabbs(s1[i]);
	s33[i] = (cmul(conj(s1[i]),s2[i])).re; 
	s43[i] = (cmul(conj(s1[i]),s2[i])).im; 
/*printf("%5.5f\t %5.5f\t %5.5f\t %5.5f\n",s11[i],s12[i],s33[i],s43[i]); 
	*/
	}

	hw 			= 7/mus; /* [cm] , maximum range in x and y for output. */
	dx 			= 2.0*hw/NN;
	dy 			= 2.0*hw/NN;

/******** MONTE CARLO *******/
 
	InitRandomGen;

/* LAUNCHNOW*/
	IT=0;/*W*/				
	QT=0;						
	UT=0;
	VT=0;
		
	IR_1=0;/*W*/				
	QR_1=0;						
	UR_1=0;
	VR_1=0;
 		
	temp=0;
		
	for (iy=0; iy<NN; iy++){
	
		for (ix=0; ix<NN; ix++) {
		
			IR[iy][ix] = 0.0;
			
			QR[iy][ix] = 0.0;
			
			UR[iy][ix] = 0.0;
			
			VR[iy][ix] = 0.0;
		}
	}	
	       
    for (jjj = 1; jjj <= 4; jjj++) {

		if (jjj == 1){
		
			S0[0] = 1;	
			S0[1] = 1;
			S0[2] = 0;
			S0[3] = 0;
			printf("launch H\n");}

		if (jjj == 2){
			
			S0[0] = 1;	
			S0[1] =	-1;
			S0[2] = 0;
			S0[3] = 0;
			printf("launch V\n");}

		if (jjj == 3){
			
			S0[0] = 1;	
			S0[1] =	0;
			S0[2] = 1;
			S0[3] = 0;
			printf("launch P\n");}

		if (jjj == 4){
			
			S0[0] = 1;	
			S0[1] =	0;
			S0[2] = 0;
			S0[3] = 1;	
			printf("launch R\n");}
       
       		
/* LAUNCH photon */
		for (i_photon = 1; i_photon <= Nphotons; i_photon++) {



/*pencil beam	*/

		x = 0.0;
		y = 0.0;
		z = 0.0;
	

/* photon direction cosines */

		U[0] = 0.0;
		U[1] = 0.0;
		U[2] = 1.0;
		

		for (i=0; i<4; i++) S[i] = S0[i]; /* set incident Stokes vector to S0 */
		for (i=0; i<4; i++) S2[i] = 0.0; /* set incident Stokes vector to S0 */
	
		photon_status = ALIVE;
		W	= 1; /* photon weight */

	
/********* ALIVE cycle *****************/
		while (photon_status == ALIVE) {
			
			
/**** HOP */
			rnd = 0; while (rnd == 0) rnd = RandomNum;  /* choose a step size */
		
			s = -log(rnd)/(mus+mua);
			x += U[0]*s;
			y += U[1]*s;
			z += U[2]*s;

	
/**** ABSORB */
		
			absorb = W*(1-albedo);
			W-= absorb;
		
	
if ( z<=0) {
		
			/*return to detector reference frame*/
	
			phi=atan2(U[1],U[0]);
			
			rotSphi(S, phi, S2);

			IR_1+=S2[0];
			QR_1+=S2[1];
			UR_1+=S2[2];
			VR_1+=S2[3];
	
			if (x >= -hw)
				ix = (int)(fabs(x + hw)/dx);
					
			if (y >= -hw)
				iy = (int)(fabs(y + hw)/dy);
					
			if (ix > MM) ix = MM; 
					
			if (iy > MM) iy = MM; 
				
			IR[iy][ix] += S2[0];
					
			QR[iy][ix] += S2[1];
					
			UR[iy][ix] += S2[2];
					
			VR[iy][ix] += S2[3]; 
				
			photon_status = DEAD;

	}
		 
	else if ( z>=slabsize) {
		
		phi=-atan2(U[1],U[0]);
		
		rotSphi(S, phi, S2);

		IT+=S2[0]*W;
		QT+=S2[1]*W;
		UT+=S2[2]*W;
		VT+=S2[3]*W;
		
		photon_status = DEAD;
						
		}/*z>slab size*/
			
	
/* SPIN */
   
/* REJECTION METHOD to choose azimuthal angle phi and deflection angle theta */
				
			do{ theta 	= acos(2*RandomNum-1);       	
			
			    phi = RandomNum*2.0*pi;              
		   							
				I0=s11[0]*S[0]+s12[0]*(S[1]*cos(2*phi)+S[2]*sin(2*phi));			
	                        
	 			ithedeg = floor(theta*nangles/pi);  
                                                 
				I=s11[ithedeg ]*S[0]+s12[ithedeg]*(S[1]*cos(2*phi)+S[2]*sin(2*phi));
				
			}while(RandomNum*I0>=I);
			
  /*------------------------------------------------------------------------------	
   Scattering : rotate to meridian plane	then scatter															
------------------------------------------------------------------------------*/
	
			updateU(U, phi, theta, U2);  /* update photon trajectory vector */                  /*sk---> important function--------------------->*/
						
			costheta=cos(theta);

			rotSphi(S, phi, S2);

			S[0]= s11[ithedeg]*S2[0]+s12[ithedeg]*S2[1];
				
			S[1]= s12[ithedeg]*S2[0]+s11[ithedeg]*S2[1];
	
			S[2]= s33[ithedeg]*S2[2]+s43[ithedeg]*S2[3];
				
			S[3]= -s43[ithedeg]*S2[2]+s33[ithedeg]*S2[3];

			temp=(sqrt(1-costheta*costheta)*sqrt(1-U2[2]*U2[2]));
			
			if ( temp==0){
				cosi=0;}
			else{
			
				if ((phi>pi) & (phi<2*pi))
					cosi=(U2[2]*costheta-U[2])/temp;	
				else
					cosi=-(U2[2]*costheta-U[2])/temp;	
				if (cosi<-1) cosi=-1;
				if (cosi>1) cosi=1;
				}
		
			sini = sqrt(1-cosi*cosi);
	
			cos22=2*cosi*cosi-1;
			
			sin22=2*sini*cosi;
			
			S2[0]=S[0];
			
			S2[1]=(S[1]*cos22-S[2]*sin22);
					
			S2[2]=(S[1]*sin22+S[2]*cos22);
					
			S2[3]=S[3];


			S[1]= S2[1]/S2[0];	
			S[2]= S2[2]/S2[0];
			S[3]= S2[3]/S2[0];
			S[0]= 1.0;
			
			for (i=0; i<3; i++) U[i] = U2[i]; /* update U */
			
/*ROULETTE*/
			rnd=0; while(rnd==0) rnd=RandomNum;
	
			if (W<THRESHOLD){
				if (rnd<=CHANCE)
					W/=CHANCE;
				else photon_status=DEAD;
			}
	
	
	} /* end of single photon launching */

}/* slab size*/

printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",IR_1/(Nphotons),QR_1/(Nphotons),UR_1/(Nphotons),VR_1/(Nphotons));	
printf("T= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",IT/(Nphotons),QT/(Nphotons),UT/(Nphotons),VT/(Nphotons));	
	
IT=0;
QT=0;
UT=0;
VT=0;
	
IR_1=0;
QR_1=0;
UR_1=0;
VR_1=0;

/* sk----> To save path of propgation */

if (jjj==1){
	target = fopen("outHI.dat","w");
	
	for (i=0; i<NN; i++) {

		fprintf(target,"%5.5f", IR[i][0]);
		
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	
	target = fopen("outHQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	target = fopen("outHU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", UR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", UR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outHV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	
	for (iy=0; iy<NN; iy++)
		for (ix=0; ix<NN; ix++) {
			IR[iy][ix] = 0.0;
			QR[iy][ix] = 0.0;
			UR[iy][ix] = 0.0;
			VR[iy][ix] = 0.0;
		}
	} 
/*111111111111111111111111111111111111111111111111111111111111111*/

	if (jjj==2) {

/* save data to file */
	target = fopen("outVI.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", IR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outVQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outVU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", UR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", UR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outVV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	for (iy=0; iy<NN; iy++)
		for (ix=0; ix<NN; ix++) {
			IR[iy][ix] = 0.0;
			QR[iy][ix] = 0.0;
			UR[iy][ix] = 0.0;
			VR[iy][ix] = 0.0;
			}
	} 
/* 222222222222222222222222222222222222222222222222222222222222222*/

	if (jjj==3) {

/* save data to file */

	target = fopen("outPI.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", IR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outPQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outPU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", UR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", UR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outPV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	for (iy=0; iy<NN; iy++)
		for (ix=0; ix<NN; ix++) {
			IR[iy][ix] = 0.0;
			QR[iy][ix] = 0.0;
			UR[iy][ix] = 0.0;
			VR[iy][ix] = 0.0;
		}
	} 
/* 33333333333333333333333333333333333333333333333333333333333333333333*/

	if (jjj==4) {
/* save data to file */
	target = fopen("outRI.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", IR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	
	target = fopen("outRQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outRU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", UR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", UR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outRV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);	

	}
	}/* end of 4 photon launchings */ 
	finish_time= clock();
	printf("Elapsed Time              = %10.2f seconds\n", (double)(finish_time-start_time)/CLOCKS_PER_SEC);

	fflush(NULL);
return 0;
} /* main routine*/ 

/*************** end MAIN ************************************/	

/*************************************************************/	
/* SUBROUTINES */
/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/


 
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double
RandomGen(char Type, long Seed, long *Status){
  static long i1, i2, ma[56];   /* ma[0] is not used. */
  long        mj, mk;
  short       i, ii;

  if (Type == 0) {              /* set seed. */
    mj = MSEED - (Seed < 0 ? -Seed : Seed);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
        mk += MBIG;
      mj = ma[ii];
    }
    for (ii = 1; ii <= 4; ii++)
      for (i = 1; i <= 55; i++) {
        ma[i] -= ma[1 + (i + 30) % 55];
        if (ma[i] < MZ)
          ma[i] += MBIG;
      }
    i1 = 0;
    i2 = 31;
  } else if (Type == 1) {       /* get a number. */
    if (++i1 == 56)
      i1 = 1;
    if (++i2 == 56)
      i2 = 1;
    mj = ma[i1] - ma[i2];
    if (mj < MZ)
      mj += MBIG;
    ma[i1] = mj;
    return (mj * FAC);
  } else if (Type == 2) {       /* get status. */
    for (i = 0; i < 55; i++)
      Status[i] = ma[i + 1];
    Status[55] = i1;
    Status[56] = i2;
  } else if (Type == 3) {       /* restore status. */
    for (i = 0; i < 55; i++)
      ma[i + 1] = Status[i];
    i1 = Status[55];
    i2 = Status[56];
  } else
    puts("Wrong parameter to RandomGen().");
  return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/************************************************************************************
 *	rotSphi(S,phi,S)
 *		Rotate S by phi [radians] and return as S
 *      multiply S for the rotational matrix of Chandrasekar or Boheren and Hoffman
 *		Uses invtan()
 ****/
void 	rotSphi(double* S, double phi, double* S2) {
	double	cos2phi, sin2phi;

	cos2phi = cos(2*phi);
	sin2phi = sin(2*phi); 
	
	S2[0] = S[0]; 
	S2[1] = S[1]*cos2phi+S[2]*sin2phi; 
	S2[2] = -S[1]*sin2phi+S[2]*cos2phi;  
	S2[3] = S[3]; 

}


/**************************************************************************
 *	updateU(U,U2)
 ****/
void 	updateU(double* U, double phi, double theta, double* U2) {
	double	ux, uy, uz, uxx, uyy, uzz, temp, sintheta, costheta, sinphi, cosphi;
	double 	pi = 3.14159265358979;
	
	ux = U[0];
	uy = U[1];
	uz = U[2];
	
	costheta = cos(theta);
	sintheta = sqrt(1.0 - costheta*costheta); 
	cosphi   = cos(phi);
	if (phi < pi)
		sinphi = sqrt(1.0 - cosphi*cosphi);   
	else
		sinphi = -sqrt(1.0 - cosphi*cosphi);
	
  /* New directional cosines. */
  if (1 - fabs(uz) <= 1.0E-12) {      /* close to perpendicular. */
    uxx = sintheta * cosphi;
    uyy = sintheta * sinphi;
    uzz = costheta * SIGN(uz);   /*  SIGN(x) is faster than division. */
    } 
  else {					/* usually use this option */
    temp = sqrt(1.0 - uz * uz);
    uxx = sintheta * (ux * uz * cosphi - uy * sinphi) / temp + ux * costheta;
    uyy = sintheta * (uy * uz * cosphi + ux * sinphi) / temp + uy * costheta;
    uzz = -sintheta * cosphi * temp + uz * costheta;
    }
  /* Update directional cosines */
  U2[0] = uxx;
  U2[1] = uyy;
  U2[2] = uzz;
}


