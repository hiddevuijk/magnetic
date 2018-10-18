#ifndef GUARD_init_random_h
#define GUARD_init_random_h


#include <vector>

template<class ran>
void init_p(std::vector<std::vector<double> >& p, ran& ranNR)
{
	int N = p.size();

	double radius = 1.;
	for(int i=0;i<N;++i) {
		double pi = acos(-1.);
		double sig = pi*ranNR.doub();
		double phi = 2*pi*ranNR.doub();
		p[i][0] = radius*sin(sig)*cos(phi);
		p[i][1] = radius*sin(sig)*sin(phi);
		p[i][2] = radius*cos(sig);

	}	
}

template<class ran>
void init_r(std::vector<std::vector<double> >& r,double L, ran& ranNR)
{

	int N = r.size();

	for(int i=0;i<N;++i){
		r[i][0] = ranNR.doub()*L;
		r[i][1] = ranNR.doub()*L;
		r[i][2] = ranNR.doub()*L;
	}
}
		


template<class ran>
void init_v(std::vector<std::vector<double> >& v, double vi, ran& ranNR)
{

	int N = v.size();

	for(int i=0;i<N;++i) {
		v[i][0] = ranNR.doub()*vi;
		v[i][1] = ranNR.doub()*vi;
		v[i][2] = ranNR.doub()*vi;
		
	}

}









#endif
