#ifndef GUARD_WALLS_H
#define GUARD_WALLS_H


#include <math.h>
#include <iostream>
#include <string>




class WallAll {
public:
	virtual XYZ f( const XYZ& r) = 0;


	// force from LJ potential
	double flj(double dist) {
		double dist6 = pow(sig/dist,6);
		return 24*eps*dist6*(2*dist6-1)/dist;
		}

	double get_sigma() const {return sig;}
	double get_epsilon() const {return eps;}
protected:
	double L;

	double sig;
	double eps;
	double rWCA;
	double L_minus_rWCA;
};

class NoWall: public WallAll {
	XYZ f(const XYZ& r)
		{ XYZ Fwall(0,0,0);
			return Fwall;
		}
	double get_sigma() const { return 0.;}
	double get_epsilon() const {return 0.;}
};

class SquareWall: public WallAll {
public:
	SquareWall() {};
	SquareWall(double sigg, double epss, double l)
		 {	sig = sigg; eps = epss;
			rWCA = sigg*pow(2.,1./6.); 
			L_minus_rWCA = l-rWCA; L = l;}

	XYZ f(const XYZ& r)
	{
		XYZ Fwall(0,0,0);
		if(r.x>L_minus_rWCA) {
			Fwall.x = -flj(L-r.x);
		} else if(r.x<rWCA) {
			Fwall.x = flj(r.x);
		}


		if(r.y>L_minus_rWCA) {
			Fwall.y = -flj(L-r.y);
		} else if(r.y<rWCA) {
			Fwall.y = flj(r.y);
		} 

		return Fwall;
	}
};


class TubeWall: public WallAll {
public:
	TubeWall() {};
	TubeWall(double sigg, double epss, double l)
		 {	sig = sigg; eps = epss;
			rWCA = sigg*pow(2.,1./6.); 
			L_minus_rWCA = l-rWCA; L = l;}

	XYZ f(const XYZ& r)
	{
		XYZ Fwall(0,0,0);
		if(r.x>L_minus_rWCA) {
			Fwall.x = -flj(L-r.x);
		} else if(r.x<rWCA) {
			Fwall.x = flj(r.x);
		}

		return Fwall;
	}
};


class Wall {
public:
	Wall();
	Wall(double sig, double eps, double L,
		std::string wallName);

	XYZ wallForce(const XYZ &r) const
		{ return wall_ptr->f(r); }

	double get_epsilon() const
		{ return wall_ptr->get_epsilon();}
	double get_sigma() const
		{return wall_ptr->get_sigma(); }
private:
	NoWall nowall;
	SquareWall squarewall;
	TubeWall tubewall;

	WallAll *wall_ptr;


};

Wall::Wall(double sig, double eps, double L,
			std::string wallName)
:	nowall(),
	squarewall(sig,eps,L),
	tubewall(sig,eps,L)
{

	if(wallName == "none") {
		wall_ptr = &nowall;
	} else if(wallName == "square") {
		wall_ptr = &squarewall;
	} else if(wallName == "tube") {
		wall_ptr = &tubewall;
	} else {
		std::cerr << "ERROR: " << wallName
			<< " is not a valid option.\n";
	}
}	

Wall::Wall()
:	nowall(),
	squarewall(0,0,0)
{
	wall_ptr = &nowall;
}

#endif
