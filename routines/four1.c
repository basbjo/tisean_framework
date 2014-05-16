#include <math.h>

#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr

void four1(data, nn, isign)
double data[];
int nn,isign;
/* See "Numerical Recipes in C" by Press, Flannery, Teukolsky, Vetterling,
 * Cambridge University Press 1988, p. 411f.
 *
 * Replaces 'data' by its discrete Fourier transform, if 'isign' is input as 1;
 * or replaces 'data' by 'nn' times its inverse discrete Fourier transform, if
 * 'isign' is input as -1. 'data' is a complex array of length 'nn', input as a
 * real array 'data[1:..2*nn]'. 'nn' MUST be an integer power of 2 (this is not
 * checked for!).
 *
 * The discrete Fourier transform of the N points h_k is
 *
 *     H_n = \sum_{k=0}^{N-1} h_k \exp\left(2\pi ikn/N\right)
 *
 * Complex array:
 * data[2*nn] = {real1, imag1, real2, imag2, ..., realnn, imagnn}
 *
 * Input array
 * +------------------+
 * | 1          real  | t = 0
 * | 2          imag  |
 * | 3          real  | t = dt
 * | 4          imag  |
 * |                  |
 * |    .             |
 * |    .             |
 * |    .             |
 * |                  |
 * | 2N - 3     real  | t = (N - 2)dt
 * | 2N - 2     imag  |
 * | 2N - 1     real  | t = (N - 1)dt
 * | 2N         imag  |
 * +------------------+
 *
 * Output array
 * +------------------+
 * | 1          real  | f = 0
 * | 2          imag  |
 * | 3          real  | f = 1/(N dt)
 * | 4          imag  |
 * |                  |
 * |    .             |
 * |    .             |
 * |                  |
 * | N - 1      real  | f = (N/2 - 1)/(N dt)
 * | N          imag  |
 * | N + 1      real  | f = +/- 1/(2 dt) (combination)
 * | N + 2      imag  |
 * | N + 2      real  | f = -(N/2 - 1)/(N dt)
 * | N + 4      imag  |
 * |                  |
 * |    .             |
 * |    .             |
 * |                  |
 * | 2N - 1      real | f = - 1/(N dt)
 * | 2N          imag |
 * +------------------+
 *
 * You can also use a routine like 'four1' without modification
 * even if your input data array is zero-offset, that is has
 * the range 'data[0..2*nn-1]'. In this case, simply decrement
 * the pointer to 'data' by one when 'four1' is invoked, e.g.
 * 'four1(data-1,1024,1;'.
 */
{
	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	//Double precision for the trigonometric reccurences
	double tempr, tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		//This is the bit-reversal section of the routine.
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >>1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	//Here begins the Danielson-Lanczos section of the routine.
	mmax = 2;
	while (n > mmax) {
		//Outer loop executed log_2 nn times.
		istep=2*mmax;
		//Initialize for the trigonometric recurrence.
		theta=6.28318530717959/(isign*mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			//Here are the two nested inner loops.
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				//This is the Danielson-Lanczos formula:
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			//Trigonometric recurrence.
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
