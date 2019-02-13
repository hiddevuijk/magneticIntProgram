#ifndef GUARD_orientation_h
#define GUARD_orientation_h

#include <vector>


void orientation(
	const std::vector<std::vector<double> >& r,
	const std::vector<std::vector<double> >& p,
	std::vector<std::vector<double> >& px,
	std::vector<double>& pxN,
	std::vector<std::vector<double> >& py,
	std::vector<double>& pyN,
	std::vector<std::vector<double> >& pz,
	std::vector<double>& pzN,
	double  bs, double L)
{

	int N = r.size();
	int jx,jy,jz;
	double x,y,z;



	for(int i=0;i<N;++i) {
		x = r[i][0] - L*floor(r[i][0]/L);
		jx = floor(x/bs);
		y = r[i][1] - L*floor(r[i][1]/L);
		jy = floor(y/bs);
		z = r[i][2] - L*floor(r[i][2]/L);
		jz = floor(z/bs);

		px[jx][0] += p[i][0];
		px[jx][1] += p[i][1];
		px[jx][2] += p[i][2];
		pxN[jx] += 1;

		py[jy][0] += p[i][0];
		py[jy][1] += p[i][1];
		py[jy][2] += p[i][2];
		pyN[jy] += 1;

		pz[jz][0] += p[i][0];
		pz[jz][1] += p[i][1];
		pz[jz][2] += p[i][2];
		pzN[jz] += 1;

	}
	
}







#endif
