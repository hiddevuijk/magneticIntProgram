#ifndef GUARD_WALLS_H
#define GUARD_WALLS_H

#include "vecmanip.h"


#include <cmath>


class Wall {
public:
	virtual void f(
		const std::vector<double>& ri,
		std::vector<double>& Fwall) = 0;

	// Lennard-Jones potential
	double ulj(double r) {
		double r6 = pow(sig/r,6);
		return 4*eps*r6*(r6-1);
		}

	// force from LJ potential
	double flj(double dist) {
		double dist6 = pow(sig/dist,6);
		return 24*eps*dist6*(2*dist6-1)/dist;
		}


protected:
	double L;
	double Ri;
	double Ro;

	double sig;
	double eps;
	double rco;
};

class NoWall: public Wall {
	void f(const std::vector<double>& ri,
		std::vector<double>& Fwall) {};
};

class Square: public Wall {
public:
	Square(double sigg, double epss,double rcoo, double l)
		 {sig = sigg; eps = epss;rco=rcoo, L = l;}
	void f( const std::vector<double>& r,
		std::vector<double>& Fwall) {
			if(r[0]>(L-rco)) {
				Fwall[0] = -flj(L-r[0]);
			} else if(r[0]<rco) {
				Fwall[0] = flj(r[0]);
			} else Fwall[0] = 0;

			if(r[1]>(L-rco)) {
				Fwall[1] = -flj(L-r[1]);
			} else if(r[1]<rco) {
				Fwall[1] = flj(r[1]);
			} else Fwall[1] = 0;

			Fwall[2] = 0;
	}
};

class TubeXpot: public Wall {
public:
	TubeXpot(double sigg, double epss,double rcoo, double l)
		 {sig = sigg; eps = epss;rco=rcoo, L = l;}
	void f( const std::vector<double>& r,
		std::vector<double>& Fwall) {
			if(r[1]>(L-rco)) {
				Fwall[1] = -flj(L-r[1]);
			} else if(r[1]<rco) {
				Fwall[1] = flj(r[1]);
			} else Fwall[1] = 0;
			Fwall[1] += 2*eps*(L/2. - r[1]); 
			Fwall[0] = 0;
			Fwall[2] = 0;
	}
};
class TubeX: public Wall {
public:
	TubeX(double sigg, double epss,double rcoo, double l)
		 {sig = sigg; eps = epss;rco=rcoo, L = l;}
	void f( const std::vector<double>& r,
		std::vector<double>& Fwall) {
			if(r[1]>(L-rco)) {
				Fwall[1] = -flj(L-r[1]);
			} else if(r[1]<rco) {
				Fwall[1] = flj(r[1]);
			} else Fwall[1] = 0;
			Fwall[0] = 0;
			Fwall[2] = 0;
	}
};

class Disk: public Wall {
public:
	Disk(double sigg, double epss, double rcoo,
				double ro, double LL) 
			{sig = sigg; eps = epss; rco=rcoo;Ro = ro;L=LL;}
	void f( const std::vector<double>& r,
		std::vector<double>& Fwall) {
		Fwall[2] = 0.;
		double d = sqrt( (r[0]-L/2)*(r[0]-L/2) +
						 (r[1]-L/2)*(r[1]-L/2) );

		if(d>Ro-rco) {
			double force = flj(Ro-d);
			Fwall[0] = -force*(r[0]-L/2)/d;
			Fwall[1] = -force*(r[1]-L/2)/d;
		} else {
			Fwall[0] = 0.;
			Fwall[1] = 0.;
		}
	}


};


class Doughnut: public Wall {
public:
	Doughnut(double sigg, double epss,double rcoo, 
				double ri,double ro, double LL) 
		{sig = sigg; eps = epss;rco=rcoo, Ri = ri; 
				Ro = ro; L = LL;}
	void f( const std::vector<double> &r,
		std::vector<double>& Fwall) {
		Fwall[2] = 0.;
		double d = sqrt( (r[0]-L/2)*(r[0]-L/2) +
						 (r[1]-L/2)*(r[1]-L/2) );
		if(d>Ro-rco) {
			double force = flj(Ro-d);
			Fwall[0] = -force*(r[0]-L/2)/d;
			Fwall[1] = -force*(r[1]-L/2)/d;
		} else if(d<Ri+rco) {
			double force = flj(d-Ri);
			Fwall[0] = force*(r[0]-L/2)/d;
			Fwall[1] = force*(r[1]-L/2)/d;
		} else {
			Fwall[0] = 0.;
			Fwall[1] = 0.;
		}
	
	}
};




#endif





