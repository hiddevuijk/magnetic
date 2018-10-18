#ifndef GUARD_integrate_h
#define GUARD_integrate_h

#include <vector>

// increments the state of the system over time T,
// in steps of dt

template<class derivobj>
void integrate(
	std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	std::vector<std::vector<double> >& v,
	std::vector<std::vector<double> >& p,
	derivobj& deriv, double T, double dt)
{
	
	double t = 0;
	while( (t+dt) < T) { 
		deriv(r,dr,v,p,dt);
		t += dt;
	}

	// integrate the remaining time
	if(t<T) 
		deriv(r,dr,v,p,T-t); 
}





#endif
