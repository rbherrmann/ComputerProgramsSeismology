#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sacsubc.h"
#include "noisemodel.h"

#define SIGN(a ) ( (a) >  0  ? (1):(-1))
#define ABS(a  ) ( (a) >  0  ? (a):-(a))
#define YES 1
#define NO  0


float gasdev(int idum);
float rstart(int iseed);
float rnor(void);
void gcmdln(int argc, char **argv);
void usage(void);
void four(float *data, int nn, int isign, float *dt, float *df);
void arr_locate(float *x, int n, float xval, int *jind);
void  initialize_noise(void);
void cleanup(void);
int npow2(int n);
float getpsd(float per, float pval, float dur);

float *x;
float *xvh, *yvh;
float *xvl, *yvl;
int npsdh, npsdl;
float *z;

/* control parameters from gcmdln */
int nout;
int npts;
float dt;
float pval;
int iseed;

/* control for random number generation */
int ifirst=0;

main(int argc, char **argv)
{
	float zv;
	int i,j,k,jind;
	float p;
	float df;
	float btime;
	float per;
	float dur;
	double xpsd;
	int nerr;
	double sums,sumn;

	/* parse command line arguments */
	gcmdln(argc, argv);

	/* initialize psd models */
	initialize_noise();

	/* initialize the random number sequence */
	zv=rstart(iseed);
	/* fill complex array with random numbers */
	for(i=0,j=0;i<npts;i++){
		zv=rnor();
		z[j++] = zv ;
		z[j++] = 0.0;
	}
		
	/* get the spectrum*/
	four(z,npts, -1, &dt, &df);
	/* normalize spectrum */
	sums = 0.0;
	sumn = 0.0;
	for(j=0 ; j <= npts/2 ; j+=2){
		sumn+=1.0;
		sums += z[j]*z[j] + z[j+1]*z[j+1] ;
	}
	sums = sqrt(sums/sumn) ;
	for(j=0 ; j < 2*npts ; j++)
		z[j] /= sums;
	four(z,(int)npts, +1, &dt, &df);
	/* the spectrum now has an RMS=1 */
	/* get the spectrum*/
	four(z,npts, -1, &dt, &df);
	/* now multiply by the noise spectrum 
         * The complex array is npts long and is in memory
         * i    Freq     R I j value
         * 0      0      0    1
         * 1      df     2    3
         * N/2-1 fn-df  N-1   N
         * N/2   fn      N   N+1
         * N/2-1 fn+df  N+2  N+3
         * N-1    -df  2N-2 2N-1
         * fn = Nyquist frequency = 1./N dt
         */
	dur = npts*dt;
	for(i=0,j=0,k=2*npts-1; i <= npts/2 ; i++){
		if(i == 0)
			per = 1./df;
		else
			per = 1./(i*df);
		/* get high model */
                xpsd=getpsd(per,pval,dur);
		z[j++] *= xpsd;
		z[j++] *= xpsd;
		/* fill in negative frequencies */
		if(i > 1){
			z[k--] *= xpsd;
			z[k--] *= xpsd;
		}
		
	}
	/* create the final time series */
	four(z,(int)npts, +1, &dt, &df);
	for(i=0,j=0; i < npts; i++, j+=2)
		x[i]=z[j];


/* make the sac file */
btime=0.0;
wsac1("O.sac",x,nout,btime,dt,&nerr);

	cleanup();
}

float gasdev(int idum)
{
	float x;
	if(ifirst ==0){
		ifirst +=  1 ;
		if(idum == 0)
			idum = 12345 ;
		x = rstart(iseed) ;
	}
	return(rnor() );
}

/* r.f -- translated by f2c (version 20100827).
*/


float rnor_0_(int n__, int iseed)
{
    /* Initialized data */

    static float aa = 12.37586f;
    static int ii = 17;
    static int jj = 5;
    static float b = .4878992f;
    static float c__ = 12.67706f;
    static float c1 = .9689279f;
    static float c2 = 1.301198f;
    static float pc = .01958303f;
    static float xn = 2.776994f;
    static float v[65] = { .340945f,.4573146f,.5397793f,.6062427f,.6631691f,
	    .7136975f,.7596125f,.8020356f,.8417227f,.8792102f,.9148948f,
	    .9490791f,.9820005f,1.0138492f,1.044781f,1.0749254f,1.1043917f,
	    1.1332738f,1.161653f,1.189601f,1.2171815f,1.2444516f,1.2714635f,
	    1.298265f,1.3249008f,1.3514125f,1.3778399f,1.4042211f,1.4305929f,
	    1.4569915f,1.4834526f,1.5100121f,1.5367061f,1.5635712f,1.5906454f,
	    1.617968f,1.6455802f,1.6735255f,1.7018503f,1.7306045f,1.7598422f,
	    1.7896223f,1.8200099f,1.851077f,1.8829044f,1.915583f,1.9492166f,
	    1.9839239f,2.019843f,2.0571356f,2.095993f,2.136645f,2.1793713f,
	    2.2245175f,2.2725185f,2.3239338f,2.3795007f,2.4402218f,2.5075117f,
	    2.5834658f,2.6713916f,2.7769943f,2.7769943f,2.7769943f,2.7769943f 
	    };
    static float u[17] = { .8668672834288f,.3697986366357f,.8008968294805f,
	    .417388977468f,.8254561579836f,.9640965269077f,.4508667414265f,
	    .6451309529668f,.164545602473f,.2787901807898f,.06761531340295f,
	    .966322633082f,.01963343943798f,.02947398211399f,.1636231515294f,
	    .3976343250467f,.2631008574685f };

    /* System generated locals */
    float ret_val, r__1, r__2;


    /* Local variables */
    static int j;
    static float s, t, x, y;
    static int ia, ib, ic, id;
    static float un;
    static int iii, jjj;
    static float vni;


/*       http://gams.nist.gov/serve.cgi/Module/NMS/RNOR/11307/ */

/* ***BEGIN PROLOGUE  RNOR */
/* ***DATE WRITTEN   810915 (YYMMDD) */
/* ***REVISION DATE  870419 (YYMMDD) */
/* ***CATEGORY NO.  L6A14 */
/* ***KEYWORDS  RANDOM NUMBERS, NORMAL DEVIATES */
/* ***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS */
/*             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U. */

/* ***PURPOSE  GENERATES NORMAL RANDOM NUMBERS, WITH MEAN ZERO AND */
/*             UNIT STANDARD DEVIATION, OFTEN DENOTED N(0,1). */
/* ***DESCRIPTION */

/*       RNOR generates normal random numbers with zero mean and */
/*       unit standard deviation, often denoted N(0,1). */
/*           From the book, "Numerical Methods and Software" by */
/*                D. Kahaner, C. Moler, S. Nash */
/*                Prentice Hall, 1988 */
/*   Use */
/*       First time.... */
/*                   Z = RSTART(ISEED) */
/*                     Here ISEED is any  n o n - z e r o  int. */
/*                     This causes initialization of the program. */
/*                     RSTART returns a real (single precision) echo of ISEED. */

/*       Subsequent times... */
/*                   Z = RNOR() */
/*                     Causes the next real (single precision) random number */
/*                           to be returned as Z. */

/* ..................................................................... */
/*                 Typical usage */

/*                    REAL RSTART,RNOR,Z */
/*                    INTEGER ISEED,I */
/*                    ISEED = 305 */
/*                    Z = RSTART(ISEED) */
/*                    DO 1 I = 1,10 */
/*                       Z = RNOR() */
/*                       WRITE(*,*) Z */
/*                 1  CONTINUE */
/*                    END */


/* ***REFERENCES  MARSAGLIA & TSANG, "A FAST, EASILY IMPLEMENTED */
/*                 METHOD FOR SAMPLING FROM DECREASING OR */
/*                 SYMMETRIC UNIMODAL DENSITY FUNCTIONS", TO BE */
/*                 PUBLISHED IN SIAM J SISC 1983. */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  RNOR */

    switch(n__) {
	case 1: goto L_rstart;
	}

/*      Load data array in case user forgets to initialize. */
/*      This array is the result of calling UNI 100000 times */
/*         with seed 305. */


/* ***FIRST EXECUTABLE STATEMENT  RNOR */

/* Fast part... */


/*   Basic generator is Fibonacci */

    un = u[(0 + (0 + (ii - 1 << 2))) / 4] - u[(0 + (0 + (jj - 1 << 2))) / 4];
    if (un < 0.f) {
	un += 1.f;
    }
    u[ii - 1] = un;
/*           U(II) and UN are uniform on [0,1) */
/*           VNI is uniform on [-1,1) */
    vni = un + un - 1.f;
    --ii;
    if (ii == 0) {
	ii = 17;
    }
    --jj;
    if (jj == 0) {
	jj = 17;
    }
/*        INT(UN(II)*128) in range [0,127],  J is in range [1,64] */
    j = (int) (u[ii - 1] * 128) % 64 + 1;
/*        Pick sign as VNI is positive or negative */
    ret_val = vni * v[j];
    if (ABS(ret_val) <= v[j - 1]) {
	return ret_val;
    }

/* Slow part; AA is a*f(0) */
    x = (ABS(ret_val) - v[j - 1]) / (v[j] - v[j - 1]);
/*          Y is uniform on [0,1) */
    y = u[ii - 1] - u[jj - 1];
    if (y < 0.f) {
	y += 1.f;
    }
    u[ii - 1] = y;
    --ii;
    if (ii == 0) {
	ii = 17;
    }
    --jj;
    if (jj == 0) {
	jj = 17;
    }

    s = x + y;
    if (s > c2) {
	goto L11;
    }
    if (s <= c1) {
	return ret_val;
    }
/* Computing 2nd power */
    r__1 = b - b * x;
    if (y > c__ - aa * exp(r__1 * r__1 * -.5f)) {
	goto L11;
    }
/* Computing 2nd power */
    r__1 = v[j];
/* Computing 2nd power */
    r__2 = ret_val;
    if (exp(r__1 * r__1 * -.5f) + y * pc / v[j] <= exp(r__2 * r__2 * -.5f)) {
	return ret_val;
    }

/* Tail part; .3601016 is 1./XN */
/*       Y is uniform on [0,1) */
L22:
    y = u[ii - 1] - u[jj - 1];
    if (y <= 0.f) {
	y += 1.f;
    }
    u[ii - 1] = y;
    --ii;
    if (ii == 0) {
	ii = 17;
    }
    --jj;
    if (jj == 0) {
	jj = 17;
    }

    x = log(y) * .3601016f;
/*       Y is uniform on [0,1) */
    y = u[ii - 1] - u[jj - 1];
    if (y <= 0.f) {
	y += 1.f;
    }
    u[ii - 1] = y;
    --ii;
    if (ii == 0) {
	ii = 17;
    }
    --jj;
    if (jj == 0) {
	jj = 17;
    }
/* Computing 2nd power */
    r__1 = x;
    if (log(y) * -2.f <= r__1 * r__1) {
	goto L22;
    }
    r__1 = xn - x;
    ret_val = r__1 * SIGN ( ret_val) ;
    return ret_val;
L11:
    r__1 = b - b * x;
    ret_val = r__1 * SIGN ( ret_val) ;
    return ret_val;


/*  Fill */

L_rstart:
    if (iseed != 0) {

/*          Set up ... */
/*              Generate random bit pattern in array based on given seed */

	ii = 17;
	jj = 5;
	ia = abs(iseed) % 32707;
	ib = 1111;
	ic = 1947;
	for (iii = 1; iii <= 17; ++iii) {
	    s = 0.f;
	    t = .5f;
/*             Do for each of the bits of mantissa of word */
/*             Loop  over 64 bits, enough for all known machines */
/*                   in single precision */
	    for (jjj = 1; jjj <= 64; ++jjj) {
		id = ic - ia;
		if (id >= 0) {
		    goto L4;
		}
		id += 32707;
		s += t;
L4:
		ia = ib;
		ib = ic;
		ic = id;
/* L3: */
		t *= .5f;
	    }
/* L2: */
	    u[iii - 1] = s;
	}
    }
/*       Return floating echo of ISEED */
    ret_val = (float) (iseed);
    return ret_val;
} /* rnor_ */

float rnor(void)
{
    return rnor_0_(0, 0);
    }

float rstart(int iseed)
{
    return rnor_0_(1, iseed);
    }

void four(float data[], int nn, int isign, float *dt, float *df){
/*
	subroutine four(data,nn,isign,dt,df)
c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(data(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  data is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     data are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array data, replacing the input data.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c-----
      real*4 data(*)
*/
	int n;
	int i, j, m, mmax, iiii, istep;
	float tempr, tempi;
	float wr, wi;
	float wstpr, wstpi;
	float sinth, theta;
	float dtt, dff;

	dtt = (*dt);
	dff = (*df);

	n = 2 * nn;
	if((dtt) == 0.0) dtt = 1./(nn*(dff)) ;
	if((dff) == 0.0) dff = 1./(nn*(dtt)) ;
	if((dtt) != (nn*(dff))) dff = 1./(nn*(dtt)) ;
	*dt = dtt;
	*df = dff;
	j = 1;
	for (i=1;i<= n ; i+=2){
		if(i < j){
			tempr = data[j-1];
			tempi = data[j  ];
			data[j-1] = data[i-1];
			data[j  ]=data[i  ];
			data[i-1] = tempr;
			data[i  ] = tempi;
		}
		m = n/2;
statement3:		if(j <= m) goto statement4 ;
			j = j-m;
			m = m/2;
			if(m >= 2)goto statement3 ;
statement4:		
		j=j+m;
    	}
	mmax = 2 ;
	while(mmax < n ){
		istep= 2 *mmax;
		theta = 6.283185307/(float)(isign*mmax);
		sinth=sin(theta/2.);
		wstpr=-2.*sinth*sinth;
		wstpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1; m <= mmax ; m +=2){
				for(i = m ; i <= n ; i+=istep){
					j=i+mmax;
					tempr=wr*data[j-1]-wi*data[j  ];
					tempi=wr*data[j  ]+wi*data[j-1];
					data[j-1]=data[i-1]-tempr;
					data[j  ]=data[i  ]-tempi;
					data[i-1]=data[i-1]+tempr;
					data[i  ]=data[i  ]+tempi;
				}
				tempr = wr;
				wr = wr*wstpr-wi*wstpi + wr;
				wi = wi*wstpr+tempr*wstpi + wi;
		}
		mmax = istep;
	}
 /*
  * 	get the correct dimensions for a Fourier Transform 
  * 	from the Discrete Fourier Transform
*/ 
	if(isign > 0){
		/*
		frequency to time domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*df);
		}
	} else {
		/*
		time to frequency domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*dt);
		}
 }
}
void arr_locate(float *x, int n, float xval, int *jind)
{
	int increase;
	int jlow, jup, jmid;
/*
 	Written by RB Herrmann Saint Louis University 2003
 	Note: the array must be ordered, however it can be ordered
 	increasing or decreasing
 
 	The technique essentially uses an interval halving technique 
 	which means that the number of comparisons is on the order of log_2(n)
 	instead of n for a linear search
 
 	Given an ordered array x[], find jind such that
 	xval belongs to ( x[jind] x[jind+1] )
 	legitimate values are 0 <= jind <= n-2
 	jind = -1 and jind = n-1 indicate failure
*/
/* do the arrays increase or decrease ? */
	
	if(x[n-1] > x[0])
		increase = YES;
	else
		increase = NO;
	/* do end member  test */
	if(increase){
		if(xval < x[0]){
			*jind = -1 ;
			return ;
		}
		if(xval > x[n-1]){
			*jind = n ;
			return ;
		}
	} else {
		if(xval > x[0]){
			*jind = 0 ;
			return ;
		}
		if(xval < x[n-1]){
			*jind = n ;
			return ;
		}
	}
	/* safety */
	if(xval == x[0]){
		*jind = 0 ;
		return ;
	}
	if(xval == x[n-1]){
		*jind = n-1 ;
		return ;
	}
	/* 
	jlow and jup are the current extremal bounds
	jmid is the current test
	*/
	jlow = 0 ;
	jup = n -1 ;
	while( ABS(jup-jlow) != 1){
		jmid = ( jlow + jup)/2 ;
		if(increase){
			if(xval > x[jmid]){
				jlow = jmid ;
			} else {
				jup  = jmid ;
			}
		} else {
			if(xval > x[jmid]){
				jup = jmid ;
			} else { 
				jlow  = jmid ;
			}
		}
	}
	*jind = jlow ;
	return ;
}

void usage(void)
{
	fprintf(stderr,
	"Usage: sacnoise -pval pval -seed seed -dt dt -npts npts \n");
	fprintf(stderr,
	"  Create time series of noise based on ASL NLNM and NHNM models. The output has units of \n");
	fprintf(stderr,
	"     m/s**2   (default)\n");
	fprintf(stderr,
	"  The noise level can be adjusted between the low and high noise models with pval\n");
	fprintf(stderr,
	"     pval=1   High noise model\n");
	fprintf(stderr,
	"     pval=0.5 mid-noise model\n");
	fprintf(stderr,
	"     pval=0   Low  noise model\n");
	fprintf(stderr," -dt   dt     (default 1.0) sample interval \n");
	fprintf(stderr," -npts npts   (default 32768) length of time series \n");
	fprintf(stderr," -pval pval   (default 0.5) \n");
	fprintf(stderr," -seed seed   (default 12345) Integer random number seed \n");
	fprintf(stderr," -h           (default false) online help\n");
	exit(EXIT_SUCCESS);
}

void gcmdln(int argc, char **argv)
{
/* parse command line arguments, e.g.,
	sacnoise -p pval 
   where
*/
	char *cp;

	/* initialize */
	pval = 0.5 ;
	iseed = -12345;
        nout = 32768;
	dt = 1.0;


	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			cp = argv[1];
			cp++;
			if(strcmp("pval",cp) == 0 ){
				argv++;
				argc--;
				cp = argv[1];
				pval = atof(cp);
			} else if(strcmp("seed",cp) == 0 ){
				argv++;
				argc--;
				cp = argv[1];
				iseed = atoi(cp);
			} else if(strcmp("dt",cp) == 0 ){
				argv++;
				argc--;
				cp = argv[1];
				dt = atof(cp);
			} else if(strcmp("npts",cp) == 0 ){
				argv++;
				argc--;
				cp = argv[1];
				nout = atoi(cp);
/*
			} else if(strcmp("L",cp) == 0 || strcmp("l",cp)==0){
				argv++;
				argc--;
				listfile = argv[1];
				have_otherpsd++;
			} else if(strcmp("A",cp) == 0 || strcmp("l",cp)==0){
				argv++;
				argc--;
				alistfile = argv[1];
				have_anotherpsd++;
			} else if(strcmp("H",cp)==0){
				window_type = HANNING;
			} else if(strcmp("S",cp)==0){
				window_type = SIN10;
			} else if(strcmp("5",cp)==0){
				do_smooth_5 = YES;
			} else if(strcmp("NT",cp)==0){
				do_title = NO;
*/
			} else if(strcmp("h",cp) == 0 || strcmp("?",cp)==0){
				usage();
			}
			argv++;
		}
	}
/* define npts for FFT as a power of 2 number that is 
 * greater than or equal to nout */
	npts = npow2(nout);
}

/* initialize arrays for the two noise models */
void  initialize_noise(void)
{
/* the noisemodel.h has a list of period PSD in order
 * 		of increasing period. The lists are terminated 
 * 		by a -1 -1. So the number of entries is 1 less than
 * 		the size of the array */
	int i;
	npsdh = sizeof(nhnm)/sizeof (struct _PXY) -1;
	xvh = (float *)calloc(npsdh, sizeof(float));
	yvh = (float *)calloc(npsdh, sizeof(float));
	for(i=0; i < npsdh ; i++){
		xvh[i] = nhnm[i].per ;
		yvh[i] = nhnm[i].psd;
	}
	npsdl = sizeof(nlnm)/sizeof (struct _PXY) -1;
	xvl = (float *)calloc(npsdl, sizeof(float));
	yvl = (float *)calloc(npsdl, sizeof(float));
	for(i=0; i < npsdl ; i++){
		xvl[i] = nlnm[i].per ;
		yvl[i] = nlnm[i].psd;
	}
	/* initialize arrays for time series */
	x = (float *)calloc(npts, sizeof(float));
	z = (float *)calloc(2*npts, sizeof(float));
}

/* free the arrays allocated by the program */
void cleanup(void)
{
	free(xvh);
	free(yvh);
	free(xvl);
	free(yvl);
	free(x);
	free(z);
}

int npow2(int n)
{
	int k, nn;
	k = 1; nn = 2;
	while(nn < n){
		nn *=2;
	}
	return(nn);

}

/* get the PSD from the noise models, average according to pval
 * then return as a factor for correcting the spectrum
 */
float getpsd(float per, float pval, float dur)
{
	int jind;
	float p, xpsdh, xpsdl, xpsd;
	/* get the value from the high noise model */
	arr_locate(xvh,npsdh,per, &jind);
	if(jind < 0){
		xpsdh = yvh[0];
	} else if(jind >= 0 && jind < npsdh -1){
		p=(per - xvh[jind])/(xvh[jind+1]-xvh[jind]);
		xpsdh = (1.-p)*yvh[jind] + p*yvh[jind+1];
	} else {
		xpsdh = yvh[npsdh-1];
	}
	/* get the value from the low nosie model */
	arr_locate(xvl,npsdl,per, &jind);
	if(jind < 0){
	xpsdl = yvl[0];
	} else if(jind >= 0 && jind < npsdl -1){
		p=(per - xvl[jind])/(xvl[jind+1]-xvl[jind]);
		xpsdl = (1.-p)*yvl[jind] + p*yvl[jind+1];
	} else {
		xpsdl = yvl[npsdl-1];
	}
	xpsd = (1.-pval)*xpsdl + pval*xpsdh;

	xpsd = pow(10.0,xpsd/20.0);
	xpsd *= sqrt(dur);
	return (xpsd);
}
