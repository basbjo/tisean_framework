#include <math.h>

void psd1(double *data,int length,double dt)
/* Replaces 'data' by its power spectral density
 * new length will be length/2+1
 */
/*Author: Bjoern Bastian Last modified: May 16, 2014 */
{
	int i,ii;
	double lastval;
	double scale=1.0*dt/length;

	//Fourier transformation
	realft(data-1,  length/2, 1);

	for (i=0,ii=0;i<(length/2+1);i++,ii+=2)
	{
		if (ii==0)
		{
			data[0] = scale * (double) pow(data[0], 2);
			lastval = scale * (double) pow(data[1], 2);
		}
		if (ii>0 && ii<length)
		{
			data[i] = scale * (double) (pow(data[ii], 2) +
                            pow(data[ii+1], 2));
		}
	}
	data[length/2] = lastval;
}
