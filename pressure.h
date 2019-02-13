#ifndef GUARD_PRESSURE_H
#define GUARD_PRESSURE_H


#include <iostream>

#include "system.h"
#include "xyz.h"

#include <fstream>
#include <string> 

class Pressure {
public:
	Pressure(double bs, double L, int N);	
	

	void sample(const System &system);

	void normalize();

	// write to out stream
	void write_bins(std::ostream &out);
		
	// write to file named outname
	void writeX(const char* outname);
	void writeY(const char* outname);
	void write_bins(const char* outname);
	//void write(const char* outname);

	unsigned int get_Nsample() {return Nsample;}
private:

	double bs;	// bin size
	unsigned int Nbin;	// number of bins

	

	std::vector<std::vector<XYZ> > Fxy;
	std::vector<double> bins;

	unsigned int Nsample;
	int N; // number of particles

};

Pressure::Pressure(double bss,double max, int Nparticles)
{
	bs = bss;
	Nbin = (unsigned int) std::ceil(max/bs);

	Fxy = std::vector<std::vector<XYZ> >(Nbin,std::vector<XYZ>(Nbin,XYZ(0,0,0)));

	bins = std::vector<double>(Nbin,0.);
	for(unsigned int i=0;i<Nbin;++i)
		bins[i] = (i+0.5)*bs;
	Nsample = 0;
	N = Nparticles; 
}	



void Pressure::sample(const System &system)
{

	++Nsample;
	// if there is a wall, calculate force
	if(system.wall.get_epsilon() != 0) {
		XYZ r;
		unsigned int jx,jy;
		for(unsigned int i=0;i<system.N;++i ) {
			r = system.r[i];
			r.pbc(system.L);
			jx = std::floor(r.x/bs);
			jy = std::floor(r.y/bs);

			if((jx<Nbin) && (jy<Nbin) ) {
				Fxy[jx][jy] += system.Fwall[i];
			}
		}
	}	
}




void Pressure::normalize() 
{
	double norm = 1./(N*Nsample);
	for(unsigned int jx = 0;jx < Nbin; ++jx ) {
		for(unsigned int jy = 0;jy < Nbin; ++jy ) {
			Fxy[jx][jy] *= norm;
		}
	}
}


void Pressure::writeX(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << Fxy[jx][jy].x;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}

void Pressure::writeY(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << Fxy[jx][jy].y;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}

void Pressure::write_bins(std::ostream &out)
{

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

}

void Pressure::write_bins(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

	out.close();

}
#endif

