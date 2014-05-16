#include <math.h>
#include <stdio.h>

void realft(data, n, isign)
double data[];
int n,isign;
/* See "Numerical Recipes in C" by Press, Flannery, Teukolsky, Vetterling,
 * Cambridge University Press 1988, p. 417f.
 *
 * Calculates the Fourier Transform of a set of '2n' real-valued data points.
 * Replaces this data (which is stored in array 'data[1..2n]') by the positive
 * frequency half of its complex Fourier Transform. The real-valued first and
 * last components of the complex transform are returned as elements 'data[1]'
 * and 'data[2]' respectively. 'n' must be a power of 2. This routine also
 * calculates the inverse transform of a complex data array if it is the
 * transform or real data. (Result in this case must be multiplied by '1/n'.)
 */
{
	int i,i1,i2,i3,i4,n2p3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	//Double precision for the trigonometric recurrences.
	double wr,wi,wpr,wpi,wtemp,theta;
	void four1();

	//Initialize the recurrence.
	theta=3.141592653589793/(double) n;
	if (isign == 1) {
		c2 = -0.5;
		//The forward transform is here.
		four1(data,n,1);
	} else {
		//Otherwise set up for an inverse transform.
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	n2p3=2*n+3;
	//Case i=1 done separately below.
	for (i=2;i<=n/2;i++) {
		i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
		//The two separate transforms are separated out of data.
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		//Here they are recombined to form the true transform of the original real data.
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		//The recurrence.
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		//Squeeze the first and last data together to get them all within the original array.
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		//This is the inverse transform for the case isign=-1.
		four1(data,n,-1);
	}
}
