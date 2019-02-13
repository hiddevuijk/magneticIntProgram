#ifndef GUARD_BFIELD_H
#define GUARD_BFIELD_H

#include <vector>


class Bfield {
public:
	virtual double f(const std::vector<double>& r,double t) = 0; 

protected:
	double B;
	double B0;
	double w;
	double wt;
	double L;
};

class BNone: public Bfield {
	double f(const std::vector<double>& r,double t)
			{return 0;}
};

class Bsigmoid: public Bfield {
public:
	Bsigmoid(double b0, double b, double l) {
			B = b; B0 = b0; L=l;}

	double f(const std::vector<double>& r, double t) {
			return B0*tanh(B*(r[1]-L/2.));}

};


class Bsaw: public Bfield {
public:
	Bsaw(double b0, double b, double l) {
			B = b; B0 = b0; L=l;}

	double f(const std::vector<double>& r, double t) {
			if(r[1]<L/4.) {
				return B0+B*0.25*L-B*(r[1]-0.25*L);	
			} else if(r[1]>3*L/4.) {
				return B0+B*0.75*L-B*(r[1]-0.75*L);
			} else {
				return B0+B*r[1];
			}
	}

};



class BsinYt: public Bfield {
public:
	BsinYt(double BB, double ww, double wtt) {
		B = BB; w=ww; wt=wtt;}

	double f(const std::vector<double>& r,double t) {
			return B*std::sin(w*(r[1]-wt*t)); }

};




class BsinY: public Bfield {
public:
	BsinY(double BB, double ww) {
		B = BB; w=ww; }

	double f(const std::vector<double>& r,double t) {
			return B*std::sin(w*r[1]); }

};

class BlinearY: public Bfield {
public:
	BlinearY(double b0, double b) {
			B = b; B0 = b0;}

	double f(const std::vector<double>& r, double t) {
			return B0+B*r[1];}

};

class BlinearR: public Bfield {
public:
	BlinearR(double b0,double b,double l) {
		B = b; B0 = b0; L = l;}

	double f(const std::vector<double>& r,double t) {
		double d = sqrt( (r[0]-L/2)*(r[0]-L/2) + 
						 (r[1]-L/2)*(r[1]-L/2));
		return B0+B*d;
	}

};







#endif

