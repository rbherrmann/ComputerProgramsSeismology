#include "sacsubc.h"
#include <stdio.h>

/* define the maximum number of points  and the float array for 
   storing the time series. Note I use calloc() in this C library
   so that array is defined in the brsac. */
#define NPTS 100
float *x;


main()
{
	int npts, nerr;
	float dt,depmax, depmin, b, e;

	int i;

	brsac(NPTS,"filt.sac", &x, &nerr);
	/* now lets get some header values */
	getfhv("DELTA",&dt,&nerr);
		printf("DELTA  :   %f\n",dt);
	getfhv("DEPMAX",&depmax,&nerr);
		printf("DEPMAX :   %f\n",depmax);
	getfhv("DEPMIN",&depmin,&nerr);
		printf("DEPMIN :   %f\n",depmin);
	/* this is necessary since the actual number of points
           may be less than the array dimension. brsac NEVER
           reads in more than the maximum allowable numebr of
           points and resets the internal header value so that
           the npts returned here never is greater than NPTS */
	getnhv("NPTS", &npts, &nerr);
		printf("NPTS   :   %f\n",npts);
        /* output the time series */
        for(i=0;i < npts ; i ++)
		printf("x[%d] = %f\n",i,x[i]);
}

