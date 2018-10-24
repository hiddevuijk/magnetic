#ifndef GUARD_BFIELD_H
#define GUARD_BFIELD_H

#include <vector>


class Bfield {
public:
	virtual double f(const std::vector<double>& r) = 0; 

protected:
	double B;
	double B0;
	double w;
	double L;
};

class BsineY: public Bfield {
public:
	BsineY(double BB, double ww) {
		B = BB; w=ww; }

	double f(const std::vector<double>& r) {
			return B*std::sin(w*r[1]); }

};

class BlinearY: public Bfield {
public:
	BlinearY(double b0, double b) {
			B = b; B0 = b0;}

	double f(const std::vector<double>& r) {
			return B0+B*r[1];}

};

class BlinearR: public Bfield {
public:
	BlinearR(double b0,double b,double l) {
		B = b; B0 = b0; L = l;}

	double f(const std::vector<double>& r) {
		double d = sqrt( (r[0]-L/2)*(r[0]-L/2) + 
						 (r[1]-L/2)*(r[1]-L/2) + 
						 (r[2]-L/2)*(r[2]-L/2) );
		return B0+B*d;
	}

};







#endif

