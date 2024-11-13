/* srotate96
	Changes:
	19 JAN 2021 - created
		Thus purpose of this code is to rotate Uz Ur Ut waveforms
		as well as strains and stresses in a clindrical coordinate system
		to a carteaian system.  The transformation matrix is

                | u_x |   | cos theta   - sin theta   0 | | u_r |
                | u_y | = | sin theta     cos theta   0 | | u_f |
                | u_z |   |    0              0       1 | | u_z |

		These are right handed coordinate systems.  Often x is north,
		y is east, and z is down. From seismological point of view, the signal
		is propagating in a direction that makes an angle of theta with respect to the
		x-axis. Thus the backazimuth is 180 + theta degrees.

		Because this is tailored to the output of strainpulse96 it will not
		permit a general set of rotations. The input is always with specific
		cylindrical component names, and the output with specific cartesian names.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sacsubc.h"
#include <ctype.h>
#include <libgen.h>

/* define function prototypes */
void gcmdln(int argc, char **argv);
void usage(char *);
void dostress_exec(void);
void dostrain_exec(void);
void dou_exec (void);
void sac_get(char *grn, float **data, int *npts) ;
void sac_put(char *grn, float *data, int npts, char* kcmpnm, float cmpinc, float cmpaz) ;

/* define useful macros and defines */
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (b) < (a) ? (a):(b) )
#define ABS(a  ) ( (a) >  0  ? (a):-(a))
#define SIGN(a ) ( (a) >  0  ? (1):(-1))
#define YES 1
#define NO  0
#define DEGRAD 0.01745329251994329576
#define LN 10000

/* global parameters */
char *protofile;
int dou, dostrain, dostress;
double az;
char mesg[100];

/* file extensions in the cyliindrical coordinate system, */
/* file extensions in the cartesian coordinate system 
	Note to avoid confusion with those above, the file 
	extension will include an underscore
	*/
/*
char* URcmp[] = { ".Ur" , ".Ut" , ".Uz" };
char* ERcmp[] = { ".Err", ".Erf", ".Erz", ".Eff", ".Efz", ".Ezz" };
char* SRcmp[] = { ".Srr", ".Srf", ".Srz", ".Sff", ".Sfz", ".Szz" };

char* URcmp[] = { "_Ux" , "_Uy" , "_Uz" };
char* ERcmp[] = { "_Exx", "_Exy", "_Exz", "_Eyy", "_Eyz", "_Ezz" };
char* SRcmp[] = { "_Sxx", "_Sxy", "_Sxz", "_Syy", "_Syz", "_Szz" };
*/

/* define the transformation matrix */
float T[3][3];


char *tstr;


int main(int argc, char **argv)
{

	int ls;
	int i,j;
	/* parse command line arguments */
	gcmdln(argc, argv);
	fprintf(stderr,"az            :%f\n",az);
        fprintf(stderr,"fileproto     :%s\n",protofile);
        fprintf(stderr,"displacement  :%d\n",dou);
        fprintf(stderr,"stress        :%d\n",dostress      );
        fprintf(stderr,"strain        :%d\n",dostrain      );

	/* define the transformation matrix */

        T[0][0] =   cos(DEGRAD*az);
	T[0][1] = - sin(DEGRAD*az) ;
	T[0][2] = 0.0 ;
	T[1][0] = - T[0][1];
        T[1][1] = T[0][0];
	T[1][2] = 0.0 ;
	T[2][0] = 0.0 ;
	T[2][1] = 0.0 ;
	T[2][2] = 1.0 ;

	/* T matrix */
	/*
	for(i=0;i<3 ; i++)
		for(j=0;j<3 ; j++)
			fprintf(stderr,"T[%d][%d] = %f\n",i,j,T[i][j]);
	*/

	/* allocate space for the input/outfile name. 
            The length is 4 more than the prototype file plus the null character */
	ls = strlen(protofile);
	tstr = calloc(ls+5,sizeof(char));
	
	
	if(dostress == YES)
		dostress_exec();
	if(dostrain == YES)
		dostrain_exec();
	if(dou == YES     )
		dou_exec();

	/* clean up */
	free(tstr);
	return 0;
}
 


void gcmdln(int argc, char **argv)
{
/* parse command line arguments, e.g.,
        srotate96 -AZ az -U -STRESS -STRAIN -FILE prototype
*/
	char *cp;

	int have_proto=0; 
	int have_az=0;
	int have_cmp=0;

	while(argc-- > 1 ) {
		if(*argv[1] == '-'){
			cp = argv[1];
			cp++;
			if(strcmp("FILE",cp) == 0 ){
				argv++;
				argc--;
				protofile = argv[1];
				have_proto++;
			} else if(strcmp("AZ",cp) == 0 ){
				argv++;
				argc--;
				az = strtod(argv[1], (char **)NULL);
				have_az++;
			} else if(strcmp("U",cp) == 0 ){
                                dou = YES;
				have_cmp++;
			} else if(strcmp("STRESS",cp) == 0 ){
				dostress = YES;
				have_cmp++;
			} else if(strcmp("STRAIN",cp) == 0 ){
				dostrain = YES;
				have_cmp++;
			} else if(strcmp("h",cp) == 0 || strcmp("?",cp)==0){
				usage(" ");
			}
			argv++;
		}
	}
	/* put in some syntax checks */
		if(have_proto == 0 ) usage("No prototype name");
		if(have_cmp == 0 ) usage("No output conversion specified");
		if(have_az == 0) usage("No azimuth specified");
}

void usage(char *ostr)
{
	fprintf(stderr,"%s\n",ostr);
	fprintf(stderr,"Usage: srotate96 -AZ az [-U|-STRESS|-STRAIN] -FILE prototype\n");
	fprintf(stderr,"  -AZ az (required) angle between r- and x-axes\n");
	fprintf(stderr,"  -FILE prototype (required) identifier for filename \n");
	fprintf(stderr,"       for the example below this could be ../NEW/005000_0100_0010\n");
	fprintf(stderr," -U  Rotate the Ur Ut Uz from strainspulse96 to Ux Uy Uz\n");
	fprintf(stderr,"       if they exist, e.g., ../NEW/005000_0100_0010.Ur etc\n");
	fprintf(stderr,"       to create 005000_0100_0010_Ux etc in the current directory\n");
	fprintf(stderr," -STRAIN  Rotate the Err Erf .. Ezz  from strainspulse96 to Exx Eyy ..\n");
	fprintf(stderr,"       if they exist, e.g., ../NEW/005000_0100_0010.Err etc\n");
	fprintf(stderr,"       to create 005000_0100_0010_Exx etc in the current directory\n");
	fprintf(stderr," -STRESS  Rotate the Srr Srf .. Szz  from strainspulse96 to Sxx Syy ..\n");
	fprintf(stderr,"       if they exist, e.g., ../NEW/005000_0100_0010.Srr etc\n");
	fprintf(stderr,"       to create 005000_0100_0010_Sxx etc in the current directory\n");

	fprintf(stderr," -h           (default false) online help\n");
	exit(EXIT_SUCCESS);
}


void dostress_exec(void)
{
	/* read the Ur Uf Uz traces and then rotate to Ux Uy and Uz */
	/* first get the cylindrical traces */
	float *Srr, *Srf, *Srz, *Sff, *Sfz, *Szz ;	/* input */
	float *Oxx, *Oxy, *Oxz, *Oyy, *Oyz, *Ozz ; /* output Cartesian */
	float x[3][3], y[3][3];

	int i,j,k,l,n,npts;

	/* read the cylindrical traces */
	sac_get(".Srr",&Srr, &npts) ;
	sac_get(".Srf",&Srf, &npts) ;
	sac_get(".Srz",&Srz, &npts) ;
	sac_get(".Sff",&Sff, &npts) ;
	sac_get(".Sfz",&Sfz, &npts) ;
	sac_get(".Szz",&Szz, &npts) ;


	Oxx = (float *)calloc(npts,sizeof(float));
	Oxy = (float *)calloc(npts,sizeof(float));
	Oxz = (float *)calloc(npts,sizeof(float));
	Oyy = (float *)calloc(npts,sizeof(float));
	Oyz = (float *)calloc(npts,sizeof(float));
	Ozz = (float *)calloc(npts,sizeof(float));

	/* rotate
	           T
	E' = T E T^
	                  T
	E'ij = Tik Xkl (T^ )lj
	     = Tik Xkl Tjl
	*/
        
	for(n=0;n < npts;n++){
		x[0][0] = Srr[n];
		x[0][1] = Srf[n];
		x[0][2] = Srz[n];
		x[1][0] = Srf[n];
		x[1][1] = Sff[n];
		x[1][2] = Sfz[n];
		x[2][0] = Srz[n];
		x[2][1] = Sfz[n];
		x[2][2] = Szz[n];

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				y[i][j] = 0.0;
				for(k=0;k<3;k++){
					for(l=0;l<3;l++){
						y[i][j] += T[i][k]*x[k][l]*T[j][l];
					}
				}
			}
		}
		Oxx[n] = y[0][0];
		Oxy[n] = y[0][1];
		Oxz[n] = y[0][2];
		Oyy[n] = y[1][1];
		Oyz[n] = y[1][2];
		Ozz[n] = y[2][2];
	}
	/* output the files, carefully defining the component and orientation */
	
	sac_put("_Sxx", Oxx, npts,"Sxx      ",90.0, 0.0) ;
	sac_put("_Sxy", Oxy, npts,"Sxy      ",90.0, 0.0) ;
	sac_put("_Sxz", Oxz, npts,"Sxz      ",90.0, 0.0) ;
	sac_put("_Syy", Oyy, npts,"Syy      ",90.0, 0.0) ;
	sac_put("_Syz", Oyz, npts,"Syz      ",90.0, 0.0) ;
	sac_put("_Szz", Ozz, npts,"Szz      ",90.0, 0.0) ;
	
	/* clean up */
	free(Srr);
	free(Srf);
	free(Srz);
	free(Sff);
	free(Sfz);
	free(Szz);

	free(Oxx);
	free(Oxy);
	free(Oxz);
	free(Oyy);
	free(Oyz);
	free(Ozz);
}

void dostrain_exec(void)
{
	/* read the Ur Uf Uz traces and then rotate to Ux Uy and Uz */
	/* first get the cylindrical traces */
	float *Err, *Erf, *Erz, *Eff, *Efz, *Ezz ;	/* input */
	float *Oxx, *Oxy, *Oxz, *Oyy, *Oyz, *Ozz ; /* output Cartesian */
	float x[3][3], y[3][3];

	int i,j,k,l,n,npts;

	/* read the cylindrical traces */
	sac_get(".Err",&Err, &npts) ;
	sac_get(".Erf",&Erf, &npts) ;
	sac_get(".Erz",&Erz, &npts) ;
	sac_get(".Eff",&Eff, &npts) ;
	sac_get(".Efz",&Efz, &npts) ;
	sac_get(".Ezz",&Ezz, &npts) ;


	Oxx = (float *)calloc(npts,sizeof(float));
	Oxy = (float *)calloc(npts,sizeof(float));
	Oxz = (float *)calloc(npts,sizeof(float));
	Oyy = (float *)calloc(npts,sizeof(float));
	Oyz = (float *)calloc(npts,sizeof(float));
	Ozz = (float *)calloc(npts,sizeof(float));

	/* rotate
	           T
	E' = T E T^
	                  T
	E'ij = Tik Xkl (T^ )lj
	     = Tik Xkl Tjl
	*/
        
	for(n=0;n < npts;n++){
		x[0][0] = Err[n];
		x[0][1] = Erf[n];
		x[0][2] = Erz[n];
		x[1][0] = Erf[n];
		x[1][1] = Eff[n];
		x[1][2] = Efz[n];
		x[2][0] = Erz[n];
		x[2][1] = Efz[n];
		x[2][2] = Ezz[n];

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				y[i][j] = 0.0;
				for(k=0;k<3;k++){
					for(l=0;l<3;l++){
						y[i][j] += T[i][k]*x[k][l]*T[j][l];
					}
				}
			}
		}
		Oxx[n] = y[0][0];
		Oxy[n] = y[0][1];
		Oxz[n] = y[0][2];
		Oyy[n] = y[1][1];
		Oyz[n] = y[1][2];
		Ozz[n] = y[2][2];
	}
	/* output the files, carefully defining the component and orientation */
	
	sac_put("_Exx", Oxx, npts,"Exx      ",90.0, 0.0) ;
	sac_put("_Exy", Oxy, npts,"Exy      ",90.0, 0.0) ;
	sac_put("_Exz", Oxz, npts,"Exz      ",90.0, 0.0) ;
	sac_put("_Eyy", Oyy, npts,"Eyy      ",90.0, 0.0) ;
	sac_put("_Eyz", Oyz, npts,"Eyz      ",90.0, 0.0) ;
	sac_put("_Ezz", Ozz, npts,"Ezz      ",90.0, 0.0) ;
	
	/* clean up */
	free(Err);
	free(Erf);
	free(Erz);
	free(Eff);
	free(Efz);
	free(Ezz);

	free(Oxx);
	free(Oxy);
	free(Oxz);
	free(Oyy);
	free(Oyz);
	free(Ozz);
}

void dou_exec (void)
{
	/* read the Ur Uf Uz traces and then rotate to Ux Uy and Uz */
	/* first get the cylindrical traces */
	float *Ur, *Uf, *Uz ;	/* input */
	float *X, *Y, *Z    ; /* output Cartesian */
	float x[3], y[3];

	int i,j,n,npts;

	/* read the cylindrical traces */
	sac_get(".Ur",&Ur, &npts) ;
	sac_get(".Ut",&Uf, &npts) ;
	sac_get(".Uz",&Uz, &npts) ;


	X = (float *)calloc(npts,sizeof(float));
	Y = (float *)calloc(npts,sizeof(float));
	Z = (float *)calloc(npts,sizeof(float));

	/* process to create the x y z traces */
	for(n=0;n < npts;n++){
		x[0] = Ur[n];
		x[1] = Uf[n];
		x[2] = Uz[n];
		for(i=0;i<3;i++){
			y[i] = 0.0;
			for(j=0;j<3;j++){
				y[i] +=   T[i][j]*x[j] ;
			}
		}
		X[n] = y[0];
		Y[n] = y[1];
		Z[n] = y[2];
	}
	/* output the files, carefully defining the component and orientation */
	
	sac_put("_Ux", X, npts,"Ux      ",90.0, 0.0) ;
	sac_put("_Uy", Y, npts,"Uy      ",90.0,90.0) ;
	sac_put("_Uz", Z, npts,"Uz      ", 0.0, 0.0) ;
	
	/* clean up */
	free(Ur);
	free(Uf);
	free(Uz);

	free(X);
	free(Y);
	free(Z);
}

void sac_get(char *grn, float **data, int *npts) 
{
	int iret;
	strcpy(tstr,protofile);
	strcat(tstr,grn);
	fprintf(stderr,"reading %s\n",tstr);
	brsac(LN, tstr, data, &iret);
	if(iret < 0){
		if(iret == -1 )
			sprintf(mesg,"%s does not exist",tstr);
		else if(iret == -2)
			sprintf(mesg,"data points exceed limit in %s ",tstr);
		else if(iret == -3)
			sprintf(mesg,"%s is not a binary SAC file ",tstr);
		else
			sprintf(mesg,"should not reach here");

		usage(mesg);
	
	}
	if(iret == -1)usage(mesg);
	getnhv("NPTS     ",npts, &iret);
}


void sac_put(char *grn, float *data, int npts, char* kcmpnm, float cmpinc, float cmpaz) 
{
	float depmin,depmax, depmen;
	float timmin, timmax;
	float dt, btime;
	int indmin, indmax;
	int nerr;
	getfhv("DELTA   ",&dt,&nerr);
	getfhv("B       ",&btime,&nerr);
	scmxmn(data,npts,&depmax,&depmin,&depmen,&indmax,&indmin) ;
        printf("Max Min : %g %g\n",depmin,depmax );
	setfhv("DEPMAX  ",depmax, &nerr) ;
	setfhv("DEPMIN  ",depmin, &nerr) ;
	setfhv("DEPMEN  ",depmen, &nerr) ;
	setfhv("TIMMAX  ",btime+indmax*dt, &nerr) ;
	setfhv("TIMMIN  ",btime+indmin*dt, &nerr) ;
	setkhv("KCMPNM  ",kcmpnm, &nerr);
	setfhv("CMPINC  ",cmpinc, &nerr);
	setfhv("CMPAZ   ",cmpaz, &nerr);
	

	strcpy(tstr,basename(protofile));
	strcat(tstr,grn);
	fprintf(stderr,"writing %s\n",tstr);
	bwsac(LN, tstr, data );
}

