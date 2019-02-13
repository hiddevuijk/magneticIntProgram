#ifndef GUARD_ORIENTATION_H
#define GUARD_ORIENTATION_H

#include "xyz.h"

class Orientation_xy {
public:
	Orientation_xy(double bs, double max);

	void sample(const System &system);

	void normalize();
 
	// write to out stream
	void writeX(std::ostream &out);
	void writeY(std::ostream &out);
	void writeZ(std::ostream &out);
	void write_bins(std::ostream &out);

	// write to file named outname
	void writeX(const char* outname);
	void writeY(const char* outname);
	void writeZ(const char* outname);
	void write_bins(const char* outname);
	void write(const char* outname);

	unsigned int get_Nsample() {return Nsample;}

private:
	double bs;
	unsigned int Nbin;

	std::vector<std::vector<XYZ> > p;
	std::vector<std::vector<unsigned int> > p_in_bin;
	std::vector<double> bins;

	unsigned int Nsample;
};




Orientation_xy::Orientation_xy(double bss,double max)
{
	bs = bss;
	Nbin = (unsigned int) std::ceil(max/bs);

	p = std::vector<std::vector<XYZ> >(Nbin,
			std::vector<XYZ>(Nbin));
	p_in_bin = std::vector<std::vector<unsigned int> >(Nbin,
			std::vector<unsigned int>(Nbin,0));

	bins = std::vector<double>(Nbin,0.);
	for(unsigned int i=0;i<Nbin;++i)
		bins[i] = (i+0.5)*bs;
	Nsample = 0;
}	



void Orientation_xy::sample(const System &system)
{
	++Nsample;
	//double x,y;
	XYZ r;
	unsigned int jx,jy;
	for(unsigned int i=0;i<system.N;++i ) {
		r = system.r[i];
		r.pbc(system.L);
		jx = std::floor(r.x/bs);
		jy = std::floor(r.y/bs);

		//x = system.r[i].x;
		//x -= system.L*std::floor(x/system.L);
		//jx = std::floor(x/bs);

		//y = system.r[i].y;
		//y -= system.L*std::floor(y/system.L);	
		//jy = std::floor(y/bs);

		if((jx<Nbin) && (jy<Nbin) ) 
			p[jx][jy] += system.p[i];
			p_in_bin[jx][jy] += 1;
	}
		
}

void Orientation_xy::normalize() 
{
	for(unsigned int jx = 0;jx < Nbin; ++jx ) {
		for(unsigned int jy = 0;jy < Nbin; ++jy ) {
			if(p_in_bin[jx][jy] > 0 ) 
				p[jx][jy] /= p_in_bin[jx][jy];
		}
	}
}

void Orientation_xy::writeX(std::ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << p[jx][jy].x;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}

void Orientation_xy::writeY(std::ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << p[jx][jy].y;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}


void Orientation_xy::writeZ(std::ostream &out)
{
	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << p[jx][jy].z;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}
}


void Orientation_xy::writeX(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << p[jx][jy].x;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}

void Orientation_xy::writeY(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << p[jx][jy].y;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}
void Orientation_xy::writeZ(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << p[jx][jy].z;
			if(jx<(Nbin-1)) out << ' '; 
		}
		if(jy<(Nbin-1)) out << '\n';
	}

	out.close();

}
void Orientation_xy::write_bins(std::ostream &out)
{

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

}

void Orientation_xy::write_bins(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int j=0;j<Nbin;++j) {
		out << bins[j];
		if(j<(Nbin-1)) out << ' '; 
	}

	out.close();

}

void Orientation_xy::write(const char* outname)
{
	std::ofstream out;
	out.open(outname);

	for(unsigned int jy=0;jy<Nbin;++jy) {
		for(unsigned int jx=0;jx<Nbin;++jx) {
			out << bins[jx] << '\t'
				<< bins[jy] << '\t'
				<< p[jx][jy].x << '\t'
				<< p[jx][jy].y << '\t'
				<< p[jx][jy].z << '\n';
		}
	}

	out.close();

}

#endif
