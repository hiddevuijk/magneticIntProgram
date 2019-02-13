#ifndef GUARD_box_muller_h
#define GUARD_box_muller_h

/*
	RandomGenerator is an object 
	with a mehod doub() that returns
	a random bouble. (nr recipes)
*/

#include <math.h>

template <typename RandomGenerator>
double ndist(RandomGenerator& ran)
{
	double x1=0.;
	double x2=0.;
	double y1=0.;
	static double y2=0.;
	static bool useY2 = false;
	double w=0.;

	// If useY2 == true, use the y2 that was 
	// generated in the previous call.

	if(useY2) {
		y1 = y2;
		useY2 = false;
	} else { 
		do {
			x1 = 2.*ran.doub() - 1.;
			x2 = 2.*ran.doub() - 1.;
			w = x1*x1 + x2*x2;
		} while(w>=1. || w == 0.);
		w = sqrt((-2.*log(w))/w);
		y1 = x1*w;
		y2 = x2*w;
		useY2 = true;
	}
	return y1;
}




#endif



