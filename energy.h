#ifndef GUARD_ENERGY_H
#define GUARD_ENERGY_H

// density of the system projected on the xy plane

#include <iostream>

#include "system.h"
#include "xyz.h"

#include <fstream>
#include <string> 

class Energy {
public:
	Energy( int Nparticles);	

	void sample(const System &system);

	void normalize();

	void write(const char* outname);

	unsigned int get_Nsample() {return Nsample;}
private:


	double Ukin;
	double Uwall;
	unsigned int Nsample;
	unsigned int N; // number of particles
};

Energy::Energy(int Nparticles)
{

	Ukin = 0.;
	Uwall = 0.;
	Nsample = 0;
	N = Nparticles;
}	



void Energy::sample(const System &system)
{

	
	for(unsigned int i=0;i<N;++i) {
		Ukin += 0.5*system.m*(system.v[i].length_sq());
		Uwall += system.wall.wallPotential(system.r[i]);
	}
	++Nsample;
		
}




void Energy::normalize() 
{

	Ukin /= Nsample*N;
	Uwall /= Nsample*N;
}

void Energy::write(const char* outname)
{
	std::ofstream out;
	out.open(outname);
	out << Ukin << ' ' << Uwall;

	out.close();

}

#endif

