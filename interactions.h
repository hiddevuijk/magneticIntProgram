#ifndef GUARD_INTERACTIONS_H
#define GUARD_INTERACTIONS_H

#include "xyz.h"

#include <math.h>

class Interactions {
public:
	Interactions() : epsilon(0.), sigma(0.), rco(0.) {}

	Interactions(double eps, double sig)
		: epsilon(eps), sigma(sig), rco(sig*pow(2.,1./6.)){}


	// calculate forces between particles in 
	// the vector r, and store in matrix F
	void get_forces(std::vector<XYZ>& F,
		const std::vector<XYZ>& r);


	void get_forces(std::vector<XYZ>& F,double& Uint,
		const std::vector<XYZ>& r,
		const std::vector<std::vector<unsigned int> >& neigh_index,
		const std::vector<unsigned int>& neigh_number);


	void get_forces(std::vector<XYZ>& F,
		const std::vector<XYZ>& r,
		const std::vector<std::vector<unsigned int> >& neigh_index,
		const std::vector<unsigned int>& neigh_number);


	double get_epsilon() const { return epsilon; };
private:

	// force between r1 and r2	
	XYZ force(const XYZ& r1,const XYZ& r2);
	// same as other, but also store potential energy
	XYZ force(const XYZ& r1,const XYZ& r2,double& U);


	double epsilon;
	double sigma;
	double rco;


	// obj. used in funcs.
	XYZ d;
	double dist, d6, f;
};


XYZ Interactions::force(const XYZ& r1,const XYZ& r2, double& Uint)
{
	Uint = 0;
	d = r1 - r2;
	dist = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	if( dist < rco) {
		d6 = sigma/dist;
		d6 = d6*d6*d6*d6*d6*d6;

		Uint = epsilon+4*epsilon*d6*(d6-1);

		f = 48*epsilon*d6*(d6-.5)/(dist*dist);
		return f*d;
	}
	return 0*d;
}	

XYZ Interactions::force(const XYZ& r1,const XYZ& r2)
{
	d = r1 - r2;
	dist = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	if( dist < rco) {
		d6 = sigma/dist;
		d6 = d6*d6*d6*d6*d6*d6;
		f = 48*epsilon*d6*(d6-.5)/(dist*dist);
		return f*d;
	}
	return 0*d;
}	


void Interactions::get_forces(
	std::vector<XYZ>& F, double& Uint, const std::vector<XYZ>& r,
	const std::vector<std::vector<unsigned int> >& neigh_index,
	const std::vector<unsigned int>& neigh_number)
{
	XYZ f;
	Uint = 0;
	double uPair;
	std::fill(F.begin(),F.end(),XYZ(0.,0.,0.));
	for(unsigned int index = 0;index<r.size();++index) {
		for(unsigned int i=0; i< neigh_number[index];++i) {
			f = force(r[index],r[neigh_index[index][i]],uPair);
			Uint += uPair;	
			F[index] += f;
			// if neigh. pair is stored once
			F[neigh_index[index][i] ] -= f;
		}
	}
}
void Interactions::get_forces(
	std::vector<XYZ>& F, const std::vector<XYZ>& r,
	const std::vector<std::vector<unsigned int> >& neigh_index,
	const std::vector<unsigned int>& neigh_number)
{
	XYZ f;
	std::fill(F.begin(),F.end(),XYZ(0.,0.,0.));
	for(unsigned int index = 0;index<r.size();++index) {
		for(unsigned int i=0; i< neigh_number[index];++i) {
			f = force(r[index],r[neigh_index[index][i]]);
			F[index] += f;
			// if neigh. pair is stored once
			F[neigh_index[index][i] ] -= f;
		}
	}
}
void Interactions::get_forces(
	std::vector<XYZ>& F, const std::vector<XYZ>& r)
{
	XYZ f;
	std::fill(F.begin(),F.end(),XYZ(0.,0.,0.));

	for(unsigned int i=0; i<(r.size()-1); ++i) {
		for(unsigned int j=i+1; j<r.size(); ++j) {
			f = force(r[i],r[j]);
			F[i] += f;
			F[j] -= f;
		}
	}
}




#endif	// GUARD_INTERACTIONS_H
