// perturbation classes

class PerturbationParameters{
public:
	enum Gauge {ORG, IRG, HG, SG}; // Outgoing Radiation Gauge (ORG), Ingoing Radiation Gauge (IRG), Harmonic Gauge (HG), and Scalar Gauge (SG)
	PerturbationParameters();
	PerturbationParameters(int s);
	PerturbationParameters(int s, Gauge gauge);
	~PerturbationParameters();
	
	int getSpinWeight();
	Gauge getGauge();
	
	void setSpinWeight(int s);
	void setGauge(Gauge gauge);
	
protected:
	int _s;
	Gauge _gauge;
};

class KerrSpacetime{
public:
	enum Coordinates {BoyerLindquist, BoyerLindquistZ};
	enum Tetrad {Carter, Kinnersley};
	KerrSpacetime();
	KerrSpacetime(double a);
	KerrSpacetime(double a, Coordinates coordinates);
	KerrSpacetime(double a, Coordinates coordinates, Tetrad tetrad);
	~KerrSpacetime();
	
	double getBlackHoleSpin();
	double getCoordinates();
	double getMetric(double radial, double polar, int i, int j);
	Complex getTetrad(double radial, double polar, int i, int j);
	
	void setBlackHoleSpin(double a);
	void setCoordinates(Coordinates coordinates);
	
protected:
	double getMetricBoyerLindquist(double r, double theta, int i, int j);
	double getMetricBoyerLindquistZ(double r, double z, int i, int j);
	double getCarterTetradBoyerLindquist(double r, double theta, int i, int j);
	double getCarterTetradBoyerLindquistZ(double r, double z, int i, int j);
	double getKinnersleyTetradBoyerLindquist(double r, double theta, int i, int j);
	double getKinnersleyTetradBoyerLindquistZ(double r, double z, int i, int j);
	double _a;
	Coordinates _coordinates;
	Tetrad _tetrad;
};

PerturbationParameters::PerturbationParameters(){
	_s = 0;
	_gauge = SG;
}
PerturbationParameters::PerturbationParameters(int s){
	_s = s;
	_gauge = SG;
}
PerturbationParameters::PerturbationParameters(int s, Gauge gauge){
	_s = s;
	_gauge = gauge;
}
PerturbationParameters::~PerturbationParameters(){}

KerrSpacetime::KerrSpacetime(){
	_a = 0;
	_coordinates = BoyerLindquist;
	_tetrad = Kinnersley;
}
KerrSpacetime::KerrSpacetime(double a){
	_a = a;
	_coordinates = BoyerLindquist;
	_tetrad = Kinnersley;
}
KerrSpacetime::KerrSpacetime(double a, Coordinates coordinates){
	_a = a;
	_coordinate = coordinates;
	_tetrad = Kinnersley;
}
KerrSpacetime::KerrSpacetime(double a, Coordinates coordinates, Tetrad tetrad){
	_a = a;
	_coordinate = coordinates;
	_tetrad = tetrad;
}
KerrSpacetime::~KerrSpacetime(){}