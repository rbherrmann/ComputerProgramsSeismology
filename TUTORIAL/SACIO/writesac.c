#include "sacsubc.h"

/* this is a sample program that creates two sac files
   one of a simple impulse and the other with the impulse passed
   through a lowpass recursive digital filter
*/

/* define the maximum number of points  and the two float arrays */
#define NPTS 100
float x[NPTS];
float y[NPTS];

void outputsac(int npts, float *arr, float dt, char *filename);

main()
{
	int i;
	float dt = 0.01 ;
	int offset = 10;
        float wc=31.415927;
        float f,a,b;

        f = 2./(dt*wc);
        a = 1. + f;
        b = 1. - f;
	/* initialize the impulse */
	for(i=0;i< NPTS;i++)
		x[i] = 0.0;
	x[offset] = 1.0/dt ;

        /* now apply a recursive digital filter to create the
           output */
        y[0] = 0.0;
        for(i=1;i < NPTS; i++){
		y[i] = (x[i]+x[i-1] - b*y[i-1])/a;
	}
	outputsac(NPTS, x, dt, "imp.sac");
	outputsac(NPTS, y, dt, "filt.sac");
}

void outputsac(int npts, float *arr, float dt, char *filename)
{
        /* create the SAC file 
           instead of using the wsac1 I will use the lower level 
           routines to provide more control on the output */
	int nerr;
	float b, e, depmax, depmin, depmen;
	/* get the extrema of the trace */
        	scmxmn(arr,npts,&depmax,&depmin,&depmen);
	/* create a new header for the new SAC file */
        	newhdr();
	/* set some header values */
        	setfhv("DEPMAX", depmax, &nerr);
        	setfhv("DEPMIN", depmin, &nerr);
        	setfhv("DEPMEN", depmen, &nerr);
        	setnhv("NPTS    ",npts,&nerr);
        	setfhv("DELTA   ",dt  ,&nerr);
                b = 0;
        	setfhv("B       ",b  ,&nerr);
        	setihv("IFTYPE  ","ITIME   ",&nerr);
        	e = b + (npts -1 )*dt;
        	setfhv("E       ",e     ,&nerr);
        	setlhv("LEVEN   ",1,&nerr);
        	setlhv("LOVROK  ",1,&nerr);
        	setlhv("LCALDA  ",1,&nerr);
	/* put is a default time for the plot */
		setnhv("NZYEAR", 1970, &nerr);
		setnhv("NZJDAY", 1, &nerr);
		setnhv("NZHOUR", 0, &nerr);
		setnhv("NZMIN" , 0, &nerr);
		setnhv("NZSEC" , 0, &nerr);
		setnhv("NZMSEC", 0, &nerr);
	/* output the SAC file */
        	bwsac(npts,filename,arr);
}
