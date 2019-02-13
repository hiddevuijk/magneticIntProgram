#ifndef GUARD_deriv_magnetic_h
#define GUARD_deriv_magnetic_h

#include <iostream>

#include <cmath>
#include <vector>
#include "vecmanip.h"
#include "box_muller.h"
#include "walls.h"
#include "bfield.h"

struct Deriv {
	public:

	Deriv(int NN, double mm, double Drr,
		double v00, Wall* wall_ptrr, 
		Bfield* bfield_ptrr, Ranq2 ranNRR):
		N(NN),m(mm),sqrt_2Dr(std::sqrt(2*Drr)), 
		v0(v00), wall_ptr(wall_ptrr),
		bfield_ptr(bfield_ptrr), ranNR(ranNRR)
		  {sqrt2 = std::sqrt(2.);}


	// evolves the state of the system to t+dt
	void  operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& v,
		std::vector<std::vector<double> >& p,
		double& t, double dt);

	private:
		int N;		// number of particles
		double m;	// mass of the particles
		double sqrt_2Dr; // sqrt(2*Dr), Dr = rot. diff. const.
		double v0;	// self-propulsion speed
		Wall* wall_ptr;
		Bfield* bfield_ptr;
		Ranq2 ranNR;	// the ranom number generator

		double sqrt2;
};



void Deriv::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& v,
		std::vector<std::vector<double> >& p,
		double& t, double dt)
{
	double sqrt_dt = std::sqrt(dt);
	// random numbers for p increment
	double etax,etay,etaz;

	double Bri; // magnetic field at position ri
	double dpx,dpy,dpz;

	// init wallForce to 0s because NoWall does noet change 
	// this vector!!!!!
	vector<double> Fwall(3,0.);
	for(int i=0;i<N;++i) {

		Bri = bfield_ptr->f(r[i],t);
		wall_ptr->f(r[i],Fwall);

		dr[i][0] = v[i][0]*dt;
		dr[i][1] = v[i][1]*dt;
		dr[i][2] = v[i][2]*dt;
		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];



		v[i][0] += (-Bri*v[i][1]*dt - v[i][0]*dt + v0*p[i][0]*dt + 
						Fwall[0]*dt + ndist(ranNR)*sqrt_dt*sqrt2)/m;
		v[i][1] += (Bri*v[i][0]*dt - v[i][1]*dt + v0*p[i][1]*dt +
						Fwall[1]*dt + ndist(ranNR)*sqrt_dt*sqrt2)/m;
		v[i][2] += (-v[i][2]*dt + v0*p[i][2]*dt + 
						Fwall[2]*dt + ndist(ranNR)*sqrt_dt*sqrt2)/m;	


		if(v0>0) { 
			etax = ndist(ranNR)*sqrt_dt*sqrt_2Dr;		
			etay = ndist(ranNR)*sqrt_dt*sqrt_2Dr;		
			etaz = ndist(ranNR)*sqrt_dt*sqrt_2Dr;		

			dpx = etay*p[i][2] - etaz*p[i][1];
			dpy = etaz*p[i][0] - etax*p[i][2];
			dpz = etax*p[i][1] - etay*p[i][0];
	
			p[i][0] += dpx;
			p[i][1] += dpy;
			p[i][2] += dpz;
			normalize(p[i]);
		}
	}
	
	t += dt;
}





#endif


