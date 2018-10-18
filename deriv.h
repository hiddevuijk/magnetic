#ifndef GUARD_deriv_magnetic_h
#define GUARD_deriv_magnetic_h

#include <iostream>

#include <cmath>
#include <vector>
#include "vecmanip.h"
#include "box_muller.h"

struct Deriv {
	public:

	Deriv(int NN, double mm, double Drr,
		double v00, double BB, double ww,
		Ranq2 ranNRR):
		N(NN),m(mm),sqrt_2Dr(std::sqrt(2*Drr)), v0(v00),
		B(BB), w(ww),ranNR(ranNRR), sqrt2(std::sqrt(2.)){}

	// evolves the state of the system to t+dt
	void  operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& v,
		std::vector<std::vector<double> >& p,
		double dt);

	private:
		int N;		// number of particles
		double m;	// mass of the particles
		double sqrt_2Dr; // sqrt(2*Dr), Dr = rot. diff. const.
		double v0;	// self-propulsion speed
		double B;	// magnetic field strength
		double w;	// w = 2 pi w0/L, w0=number of periods
		double sqrt2;
		Ranq2 ranNR;	// the ranom number generator

		double Br(const std::vector<double>& ri)
			{ return B*std::sin(w*ri[1]);}
		
		
};



void Deriv::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& v,
		std::vector<std::vector<double> >& p,
		double dt)
{
	double sqrt_dt = std::sqrt(dt);

	// random numbers for p increment
	double etax,etay,etaz;

	double Bri; // magnetic field at position ri
	double norm, norm2; // used to normalize v vector
	double vx, vy;  // values used in clalculation of v incr.
	double dpx,dpy,dpz;
	double v0x,v0y,v0z;
	for(int i=0;i<N;++i) {

		Bri = Br(r[i]);

		dr[i][0] = v[i][0]*dt;
		dr[i][1] = v[i][1]*dt;
		dr[i][2] = v[i][2]*dt;
		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];

		v0x = v[i][0];
		v0y = v[i][1];
		v0z = v[i][2];

		//norm = sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]);
		vx = v[i][0] - Bri*v[i][1]*dt/m;
		vy = v[i][1] + Bri*v[i][0]*dt/m;
		norm2 = sqrt(vx*vx + vy*vy);
		//vx *= norm/norm2;
		//vy *= norm/norm2;

		v[i][0] = vx + (-v[i][0]*dt+v0*p[i][0]*dt+ndist(ranNR)*sqrt_dt*sqrt2)/m;
		v[i][1] = vy + (-v[i][1]*dt+v0*p[i][1]*dt+ndist(ranNR)*sqrt_dt*sqrt2)/m;
		v[i][2] += (-v[i][2]*dt+v0*p[i][2]*dt+ndist(ranNR)*sqrt_dt*sqrt2)/m;	


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
}





#endif


