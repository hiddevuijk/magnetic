
class Wall {
public:
	virtual void f(
		const std::vector<double>& vi,
		std::vector<double>& Fwall);

protected:
	double L;
	double Ri;
	double Ro;

	double sig;
	double eps;

};

class NoWall: public Wall {
	void f(const std::vector<double>& vi,
		std::vector<double>& Fwall) {};
}
class TubeX: public Wall {
public:
	TubeX(double sigg, double epss,double l)
		 {sig = sigg; eps = epss; L = l;}
	void f( const std::vector<double>& vi,
		std::vector<double>& Fwall);
};

class Disk: public Wall {
public:
	Disk(double sigg, double epss, double ro) 
			{sig = sigg; eps = epss; Ro = ro;}
	void f( const std::vector<double>& vi,
		std::vector<double>& Fwall);
};


class Doughnut: public Wall {
public:
	Doughnut(double sigg, double epss, double ri,double ro) 
		{sig = sigg; eps = epss; Ri = ri; Ro = ro;}
	void f( const std::vector<double> &vi,
		std::vector<double>& Fwall};
}










