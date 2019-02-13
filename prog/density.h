#ifndef GUARD_density_h
#define GUARD_density_h

#include <vector>

// calculates the density in the xy plane
// it adds the number of particles in each bin in state r
// to rho, so the normalization should be done later
// bs is the binsize
// L is the box size
std::vector<std::vector<double> > density(
	const std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& rho,
	double bs,double L, int nbin)
{


	int N = r.size();
	int jx,jy;
	double x,y;

	for(int i=0;i<N;++i) {
		x = r[i][0] - L*floor(r[i][0]/L);
		y = r[i][1] - L*floor(r[i][1]/L);

		jx = floor(x/bs);
		jy = floor(y/bs);
		rho[jx][jy] += 1.;
	}
	return rho;
}


#endif


