#ifndef GUARD_flux_h
#define GUARD_flux_h


#include <vector>

void fluxXY(
	std::vector<std::vector<double> >& fluxX,
	std::vector<std::vector<double> >& fluxY,
	const std::vector<std::vector<double> >& r,
	const std::vector<std::vector<double> >& dr,
	const std::vector<std::vector<double> >& v,
	double bs, double L,int nbin)
{

	double x,y;
	int jx,jy;
	int N = r.size();


	// use r - a*dr as position variable
	double a = 0.;

	for(int i=0;i<N;++i) {
		x = r[i][0] - a*dr[i][0];
		x = x - L*floor(x/L);
		jx = floor(x/bs);

		y = r[i][1] - a*dr[i][0];
		y = y - L*floor(y/L);
		jy = floor(y/bs);

		fluxX[jx][jy] += v[i][0];
		fluxY[jx][jy] += v[i][1];
	}
}





#endif

