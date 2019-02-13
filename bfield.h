#ifndef GUARD_BFIELD_H
#define GUARD_BFIELD_H

#include "xyz.h"

#include <iostream>
#include <vector>

class BfieldAll {
public:
	virtual double f(const XYZ& r) = 0; 

protected:
	double B;
	double w;
	double L;
};

class BNone: public BfieldAll {
public:
	BNone() {B = 0;}
	double f(const XYZ& r)
			{return 0;}
};


class BSine: public BfieldAll {
public:
	BSine() { B = 0.; }
	BSine(double BB, double ww) {
		B = BB; w=ww; }

	double f(const XYZ& r) {
			return B*std::sin(w*r.x); }

};

class Bhill: public BfieldAll {
public:
	Bhill() { B = 0;}
	Bhill(double BB, double LL) {
		B = BB; L = LL; }

	double f(const XYZ& r) {
		if( (r.x < 0.2*L) or (r.x>0.8*L) ){
			return 0.;
		} else if( r.x<0.5*L) {
			return B*(r.x-0.2*L);
		} else {
			return B*(0.8*L-r.x);
		} 
	}
};

class Blin: public BfieldAll {
public:
	Blin() { B = 0;}
	Blin(double BB) {B = BB; }

	double f(const XYZ& r) {
		return r.x*B;
	}
};


class BlinAsym: public BfieldAll {
public:
	BlinAsym() { B = 0;}
	BlinAsym(double BB,double LL) {B = BB; L=LL;}

	double f(const XYZ& r) {
		return r.x*B - B*0.5*L;
	}
};



class Bfield {
public:
	Bfield();

	Bfield(double B, double w, double L,
			std::string BType);

	double get_field(const XYZ& r) const
		{return bfield_ptr->f(r);}

private:
	BNone bnone;
	BSine bsine;
	Bhill bhill;
	Blin  blin;
	BlinAsym blinasym;
	
	BfieldAll *bfield_ptr;
};



Bfield::Bfield()
:
	bnone(), bsine(), bhill()
{
	bfield_ptr = &bnone;
}


Bfield::Bfield(double B, double w, double L,
			std::string BType)
:
	bnone(), bsine(B,w), bhill(B,L),
	blin(B), blinasym(B,L)
{
	if( BType == "none") {
		bfield_ptr = &bnone;	
	} else if( BType == "sine" ) {
		bfield_ptr = &bsine;
	} else if( BType == "hill" ) {
		bfield_ptr = &bhill;
	} else if( BType == "linear") {
		bfield_ptr = &blin;
	} else if( BType == "linearAsym") {
		bfield_ptr = & blinasym;
	} else {
		std::cerr << "ERROR: " << BType
			<< " is not a valid option.\n";
	}


}

#endif

