// radialsolver.c

#include "radialsolver.hpp"

#define TEUK_ODE_REL_ERR 1.e-14
#define TEUK_ODE_STEP_SIZE_INIT 1.e-8
#define DELTA_R_HORIZON 1.e-2
#define R_INFINITY 50.
#define BOOST_INTEGRATE 0
#define MST_AUTO_FREQ 1.e-2
#define MST_AUTO_RADIAL 1000.
#define MST_AUTO_POINTS_MAX 256
#define ZERO_FREQ_MAX 1.e-11
#define ASYM_AUTO_RMIN 100

// Constructor
RadialTeukolsky::RadialTeukolsky(double a, int s, int L, int m, double omega, const Vector r):
	_a(a), _s(s), _L(L), _m(m), _omega(omega), _radialPoints(r), _horizonBoundary(0.), _infinityBoundary(0.),
	_horizonBoundarySolution(0., DBL_EPSILON), _horizonBoundaryDerivative(0., DBL_EPSILON), _infinityBoundarySolution(0., DBL_EPSILON), _infinityBoundaryDerivative(0., DBL_EPSILON), _inSolution(r.size(), 0.), _inDerivative(r.size(), 0.), _upSolution(r.size(), 0.), _upDerivative(r.size(), 0.)
{
	_lambda = swsh_eigenvalue(_s, _L, _m, _omega*_a);
}
RadialTeukolsky::RadialTeukolsky(double a, int s, int L, int m, double omega, double lambda, const Vector r):
	_a(a), _s(s), _L(L), _m(m), _omega(omega), _lambda(lambda), _radialPoints(r), _horizonBoundary(0.), _infinityBoundary(0.),
	_horizonBoundarySolution(0., DBL_EPSILON), _horizonBoundaryDerivative(0., DBL_EPSILON), _infinityBoundarySolution(0., DBL_EPSILON), _infinityBoundaryDerivative(0., DBL_EPSILON), _inSolution(r.size(), 0.), _inDerivative(r.size(), 0.), _upSolution(r.size(), 0.), _upDerivative(r.size(), 0.)
{if(_radialPoints.size() == 0){std::cout << "(RADIALSOLVER) RadialTeukolsky failed to initialize.\n";} }
// Destructor
RadialTeukolsky::~RadialTeukolsky() {}

// Access parameters that initialize class
double RadialTeukolsky::getBlackHoleSpin(){ return _a; }
int RadialTeukolsky::getSpinWeight(){ return _s; }
int RadialTeukolsky::getSpheroidalModeNumber(){ return _L; }
int RadialTeukolsky::getAzimuthalModeNumber(){ return _m; }
double RadialTeukolsky::getModeFrequency(){ return _omega; }
double RadialTeukolsky::getSpinWeightedSpheroidalEigenvalue(){ return _lambda; }

// Generate retarded boundary conditions for Teukolsky equation. This is neccessary if one is
// going to generate solutions with any of the integration methods {HBL, GSN, TEUK}
int RadialTeukolsky::ASYMsolveBoundary(BoundaryCondition bc){
	int success = 1;
	int boundaryFlagMax = 40;
	if( bc == In ){
		teuk_in_ASYM_series(_horizonBoundarySolution, *this, _horizonBoundary);
	 	teuk_in_derivative_ASYM_series(_horizonBoundaryDerivative, *this, _horizonBoundary);

		int boundaryFlag = 0;
		double rplus = (1. + sqrt(1. - _a*_a));
		while( (std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			_horizonBoundary = rplus + 0.5*(_horizonBoundary - rplus);
			teuk_in_ASYM_series(_horizonBoundarySolution, *this, _horizonBoundary);
			teuk_in_derivative_ASYM_series(_horizonBoundaryDerivative, *this, _horizonBoundary);
			boundaryFlag++;
		}

		if((std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			success = 0; // consider method failed
		}
	}else{
		teuk_up_ASYM_series(_infinityBoundarySolution, *this, _infinityBoundary);
		teuk_up_derivative_ASYM_series(_infinityBoundaryDerivative, *this, _infinityBoundary);

		int boundaryFlag = 0;
		while( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			_infinityBoundary *= 2;
			teuk_up_ASYM_series(_infinityBoundarySolution, *this, _infinityBoundary);
			teuk_up_derivative_ASYM_series(_infinityBoundaryDerivative, *this, _infinityBoundary);
			boundaryFlag++;
		}

		if( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			success = 0; // consider method failed
		}
	}

	return success;
}

int RadialTeukolsky::MSTsolveBoundary(BoundaryCondition bc){
	int success = 1;
	int boundaryFlagMax = 40;
	if( bc == In ){
		teuk_in_MST_series_boundary(_horizonBoundarySolution, _horizonBoundaryDerivative, *this, _horizonBoundary);

		int boundaryFlag = 0;
		double rplus = 1. + sqrt(1. - _a*_a);
		while( (std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || std::abs(_horizonBoundarySolution.getValue()) == 0. || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			_horizonBoundary = rplus + 0.5*(_horizonBoundary - rplus);
			teuk_in_MST_series_boundary(_horizonBoundarySolution, _horizonBoundaryDerivative, *this, _horizonBoundary);
			boundaryFlag++;
		}

		if((std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || std::abs(_horizonBoundarySolution.getValue()) == 0. || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			success = 0; // consider method failed
		}
	}else{
		teuk_up_MST_series_boundary(_infinityBoundarySolution, _infinityBoundaryDerivative, *this, _infinityBoundary);

		int boundaryFlag = 0;
		while( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || std::abs(_infinityBoundarySolution.getValue()) == 0. || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			_infinityBoundary *= 2;
			teuk_up_MST_series_boundary(_infinityBoundarySolution, _infinityBoundaryDerivative, *this, _infinityBoundary);
			boundaryFlag++;
		}

		if( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || std::abs(_infinityBoundarySolution.getValue()) == 0. || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && boundaryFlag < boundaryFlagMax){
			success = 0; // consider method failed
		}
	}

	return success;
}

void RadialTeukolsky::solveBoundaryPoint(BoundaryCondition bc){
	if( bc == In ){
		double rplus = 1. + sqrt(1. - _a*_a);
		_horizonBoundary = 1.1*rplus;
		if(_horizonBoundary >= _radialPoints.front()){
			_horizonBoundary = rplus + 0.5*(_radialPoints.front() - rplus);
		}
	}else{
		_infinityBoundary = R_INFINITY;
		if(20./std::abs(_omega) > _infinityBoundary && std::abs(_omega) > DBL_EPSILON){
			_infinityBoundary = 20./std::abs(_omega);
			if(_infinityBoundary > 50*_radialPoints.back()){
				_infinityBoundary = 50*_radialPoints.back();
			}
		}
		if(_infinityBoundary <= _radialPoints.back()){
			_infinityBoundary = 1.1*_radialPoints.back();
		}
	}
}

void RadialTeukolsky::generateRetardedBoundaryCondition(BoundaryCondition bc, SolutionMethod method){
	solveBoundaryPoint(bc);
	int test = 1;
	if( method != ASYM && method != MST ){
		test = ASYMsolveBoundary(bc);
		if( test == 0 && validMethodDomain(MST) ){
			test = MSTsolveBoundary(bc);
		}

		if( test == 0 ){
			ASYMsolveBoundary(bc);
		}
	}else if( method == ASYM ){
		test = ASYMsolveBoundary(bc);
	}else if( method == MST ){
		test = MSTsolveBoundary(bc);
	}

	if( test == 0 ){
		if(bc == Up){
			std::cout << "(RADIALSOLVER) ERROR: Calculation of infinity side boundary value failed for L = "<<_L<<", m = "<<_m<<", omega = "<<_omega<<". Attempted to push boundary out to r = "<<_infinityBoundary<<" \n";
		}else{
			std::cout << "(RADIALSOLVER) ERROR: Calculation of horizon side boundary value failed for L = "<<_L<<", m = "<<_m<<", omega = "<<_omega<<". Attempted to push boundary out to r = "<<_horizonBoundary<<" \n";
		}
	}
}

void RadialTeukolsky::generateRetardedBoundaryConditions(SolutionMethod method){
	generateRetardedBoundaryCondition(In, method);
	generateRetardedBoundaryCondition(Up, method);
}

// void RadialTeukolsky::generateRetardedBoundaryConditions(SolutionMethod method){
// 	_horizonBoundary = 1. + 1.1*sqrt(1. - _a*_a);
// 	_infinityBoundary = R_INFINITY;
// 	int flagHor = 0, flagInf = 0;

// 	// if MST method is specified, flag this option
// 	if(method == MST){
// 		flagHor = 1;
// 		flagInf = 1;
// 	}

// 	// if AUTO or ASYM methods selected, first calculate boundary data using series expansions
// 	// around singular points
// 	if( flagHor == 0 && flagInf == 0 ){
// 		if(20./std::abs(_omega) > _infinityBoundary && std::abs(_omega) > DBL_EPSILON){
// 			_infinityBoundary = 20./std::abs(_omega);
// 			if(_infinityBoundary/_radialPoints.back() > 50){
// 				_infinityBoundary = 50*_radialPoints.back();
// 			}
// 		}
// 		// check to see that the boundary is at least placed a small distance away from the first point we want to store
// 		if(_radialPoints.back() + 0.25 > _infinityBoundary){
// 			_infinityBoundary = _radialPoints.back() + 20;
// 		}

// 		teuk_in_ASYM_series(_horizonBoundarySolution, *this, _horizonBoundary);
// 	 	teuk_in_derivative_ASYM_series(_horizonBoundaryDerivative, *this, _horizonBoundary);
// 		teuk_up_ASYM_series(_infinityBoundarySolution, *this, _infinityBoundary);
// 		teuk_up_derivative_ASYM_series(_infinityBoundaryDerivative, *this, _infinityBoundary);

// 		// if series expansion failed, try pushing out the boundary
// 		if( std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue())) ){
// 			_infinityBoundary *= 2;
// 			teuk_up_ASYM_series(_infinityBoundarySolution, *this, _infinityBoundary);
// 			teuk_up_derivative_ASYM_series(_infinityBoundaryDerivative, *this, _infinityBoundary);
// 		}
// 		if( std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))){
// 			_horizonBoundary = 0.5*(_horizonBoundary - (1. + sqrt(1. - _a*_a)));
// 			teuk_in_ASYM_series(_horizonBoundarySolution, *this, _horizonBoundary);
// 			teuk_in_derivative_ASYM_series(_horizonBoundaryDerivative, *this, _horizonBoundary);
// 		}

// 		// if AUTO method selected, test whether these series expansions converged. If they did not,
// 		// then switch flags to use MST expansions
// 		// if different method specified, keep pushing out boundaries until series expansions converge
// 		if(method == AUTO && _omega < 0.3){
// 			int boundaryFlag = 0;
// 			// first try moving out the boundary twice to see if we can get a valid asymptotic expansion
// 			while( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && boundaryFlag < 4){
// 				_infinityBoundary *= 2;
// 				teuk_up_ASYM_series(_infinityBoundarySolution, *this, _infinityBoundary);
// 				teuk_up_derivative_ASYM_series(_infinityBoundaryDerivative, *this, _infinityBoundary);
// 				boundaryFlag++;
// 			}
// 			boundaryFlag = 0;
// 			while( (std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && boundaryFlag < 4){
// 				_horizonBoundary = 0.5*(_horizonBoundary - (1. + sqrt(1. - _a*_a)));
// 				teuk_in_ASYM_series(_horizonBoundarySolution, *this, _horizonBoundary);
// 			 	teuk_in_derivative_ASYM_series(_horizonBoundaryDerivative, *this, _horizonBoundary);
// 				boundaryFlag++;
// 			}

// 			// if that did not work, now try the MST method to generate boundary data
// 			if( std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue())) ){
// 				flagInf = 1;
// 			}
// 			if( std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue())) ){
// 				flagHor = 1;
// 			}
// 		}else{
// 			while( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && _infinityBoundary < pow(10, 6)*_radialPoints.back()){
// 				_infinityBoundary *= 2;
// 				teuk_up_ASYM_series(_infinityBoundarySolution, *this, _infinityBoundary);
// 				teuk_up_derivative_ASYM_series(_infinityBoundaryDerivative, *this, _infinityBoundary);
// 			}
// 			while( (std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && _horizonBoundary > (1. + sqrt(1. - _a*_a)) + 10.*DBL_EPSILON){
// 				_horizonBoundary = 0.5*(_horizonBoundary - (1. + sqrt(1. - _a*_a)));
// 				teuk_in_ASYM_series(_horizonBoundarySolution, *this, _horizonBoundary);
// 			 	teuk_in_derivative_ASYM_series(_horizonBoundaryDerivative, *this, _horizonBoundary);
// 			}
// 		}
// 		// std::cout << "Horizon side data from asymptotics = " << _horizonBoundarySolution << " for r = " << _horizonBoundary << "\n";
// 		// std::cout << "Infinity side data from asymptotics = " << _infinityBoundarySolution << " for r = " << _infinityBoundary << "\n";
// 		// std::cout << "Horizon side derivative data from asymptotics = " << _horizonBoundaryDerivative << " for r = " << _horizonBoundary << "\n";
// 		// std::cout << "Infinity side derivative data from asymptotics = " << _infinityBoundaryDerivative << " for r = " << _infinityBoundary << "\n";
// 	}

// 	// these flags tell us which boundary data to calculate with MST expansions
// 	if( flagHor && flagInf && std::abs(_omega) > ZERO_FREQ_MAX ){
// 		std::cout << "(RADIALSOLVER) Using MST expansions to generate boundary data. \n";
// 		if(std::abs(_omega) < 1. && std::abs(2.*_omega*_radialPoints.back()) >= 0.5){
// 			_infinityBoundary = 20./std::abs(_omega);
// 		}else if(std::abs(2.*_omega*_radialPoints.back()) < 0.5){
// 			_infinityBoundary = 0.5/std::abs(_omega);
// 		}

// 		teuk_MST_series_boundary(_horizonBoundarySolution, _horizonBoundaryDerivative, _infinityBoundarySolution, _infinityBoundaryDerivative, *this, _horizonBoundary, _infinityBoundary);

// 	}else if(flagHor && std::abs(_omega) > ZERO_FREQ_MAX){
// 		std::cout << "(RADIALSOLVER) Generating horizon side MST expansions. \n";
// 		teuk_in_MST_series_boundary(_horizonBoundarySolution, _horizonBoundaryDerivative, *this, _horizonBoundary);
// 	}else if(flagInf && std::abs(_omega) > ZERO_FREQ_MAX){
// 		std::cout << "(RADIALSOLVER) Generating infinity side MST expansions. \n";
// 		if(std::abs(_omega) < 1. && std::abs(2.*_omega*_radialPoints.back()) >= 0.5){
// 			_infinityBoundary /= std::abs(_omega);
// 		}else if(std::abs(2.*_omega*_radialPoints.back()) < 0.5){
// 			_infinityBoundary = 0.5/std::abs(_omega);

// 		}
// 		teuk_up_MST_series_boundary(_infinityBoundarySolution, _infinityBoundaryDerivative, *this, _infinityBoundary);
// 	}else{
// 		// std::cout << "(RADIALSOLVER) Using series expansions near singular points to generate boundary data. \n";
// 	}

// 	_inTransmissionAmplitude = 1.;
// 	_inIncidenceAmplitude = 0.;
// 	_inReflectionAmplitude = 0.;
// 	_upTransmissionAmplitude = 1.;
// 	_upIncidenceAmplitude = 0.;
// 	_upReflectionAmplitude = 0.;

// 	while( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && _infinityBoundary < pow(10, 6)){
// 		_infinityBoundary *= 2;
// 		teuk_up_ASYM_series(_infinityBoundarySolution, *this, _infinityBoundary);
// 		teuk_up_derivative_ASYM_series(_infinityBoundaryDerivative, *this, _infinityBoundary);
// 	}
// 	while( (std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || (std::abs(_horizonBoundarySolution.getValue())) == 0. || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && _horizonBoundary > (1. + sqrt(1. - _a*_a)) + 10.*DBL_EPSILON){
// 		_horizonBoundary = 0.5*(_horizonBoundary - (1. + sqrt(1. - _a*_a)));
// 		teuk_in_ASYM_series(_horizonBoundarySolution, *this, _horizonBoundary);
// 		teuk_in_derivative_ASYM_series(_horizonBoundaryDerivative, *this, _horizonBoundary);
// 	}

// 	if( (std::abs(_infinityBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundarySolution.getValue())) || isinf(std::abs(_infinityBoundarySolution.getValue())) || std::abs(_infinityBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_infinityBoundaryDerivative.getValue())) || isinf(std::abs(_infinityBoundaryDerivative.getValue()))) && _infinityBoundary > pow(10, 6)){
// 		std::cout << "(RADIALSOLVER) ERROR: Calculation of infinity side boundary value failed for L = "<<_L<<", m = "<<_m<<", omega = "<<_omega<<". Attempted to push boundary out to r = "<<_infinityBoundary<<" \n";
// 	}
// 	if( (std::abs(_horizonBoundarySolution.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundarySolution.getValue())) || isinf(std::abs(_horizonBoundarySolution.getValue())) || std::abs(_horizonBoundaryDerivative.getPrecision()) > 1.e-13 || isnan(std::abs(_horizonBoundaryDerivative.getValue())) || isinf(std::abs(_horizonBoundaryDerivative.getValue()))) && _horizonBoundary < (1. + sqrt(1. - _a*_a)) + 10.*DBL_EPSILON){
// 		std::cout << "(RADIALSOLVER) ERROR: Calculation of horizon side boundary value failed for L = "<<_L<<", m = "<<_m<<", omega = "<<_omega<<". Attempted to push boundary out to r = "<<_horizonBoundary<<" \n";
// 	}

// 	//std::cout << "Horizon side data = " << _horizonBoundarySolution << " for r = " << _horizonBoundary << "\n";
// 	//std::cout << "Infinity side data = " << _infinityBoundarySolution << " for r = " << _infinityBoundary << "\n";
// 	// std::cout << std::setprecision(6);
// }

// alternatively one can specify their own boundary conditions
void RadialTeukolsky::setBoundaryConditions(BoundaryCondition bc, Complex R, Complex Rp, double r){
	if(bc == In){
		_horizonBoundary = r;
		_horizonBoundarySolution = Result(R, DBL_EPSILON);
		_horizonBoundaryDerivative = Result(Rp, DBL_EPSILON);
	}else{
		_infinityBoundary = r;
		_infinityBoundarySolution = Result(R, DBL_EPSILON);
		_infinityBoundaryDerivative = Result(Rp, DBL_EPSILON);
	}
}

// Check to see if the solution method failed to generate non-trivial and finite solutions.
// Returns true if solution method failed, false if successful
bool RadialTeukolsky::failCheck(){
	for(size_t i = 0; i < _inSolution.size(); i++){
		if(isnan(std::abs(_inSolution[i])) || std::abs(_inSolution[i]) == 0. || isnan(std::abs(_upSolution[i])) || std::abs(_upSolution[i]) == 0. || isnan(std::abs(_inDerivative[i])) || isnan(std::abs(_upDerivative[i]))){
			// std::cout << "Asymptotic series solution failed for r = "<< _radialPoints[i] <<" \n";
			return true;
		}
	}
	return false;
}

// Check to see if the solution method failed to generate non-trivial and finite solutions for a given boundary condition
bool RadialTeukolsky::failCheck(BoundaryCondition bc){
	if(bc == In){
		for(size_t i = 0; i < _inSolution.size(); i++){
			if(isnan(std::abs(_inSolution[i])) || std::abs(_inSolution[i]) == 0. || isnan(std::abs(_inDerivative[i]))){
				// std::cout << "Asymptotic series solution failed for r = "<< _radialPoints[i] <<" \n";
				return true;
			}
		}
		return false;
	}else{
		for(size_t i = 0; i < _upSolution.size(); i++){
			if(isnan(std::abs(_upSolution[i])) || std::abs(_upSolution[i]) == 0. || isnan(std::abs(_upDerivative[i]))){
				// std::cout << "Asymptotic series solution failed for r = "<< _radialPoints[i] <<" \n";
				return true;
			}
		}
		return false;
	}
}


// Helper function to see if particular method is expected to converge in an area
// of parameter space in which we want to evaluate
bool RadialTeukolsky::validMethodDomain(SolutionMethod method){
	if(method == ASYM){
		return (std::abs(_omega*_radialPoints.front()) > ASYM_AUTO_RMIN);
	}else if(method == MST){
		return ((std::abs(_radialPoints.back()*_omega) < MST_AUTO_FREQ && _radialPoints.size() < MST_AUTO_POINTS_MAX) || std::abs(10.*_radialPoints.back()*_omega) < MST_AUTO_FREQ);
	}else{
		return true;
	}
}

// Helper function for solving the Teukolsky equation using hyperboloidal
// slicing. The make_stable option indicates whether or not we only 
// integrate the spin-weight solutions that are numerically stable for a given
// boundary condition
void RadialTeukolsky::HBLsolve(BoundaryCondition bc, bool make_stable){
	if(bc == In){
		if(_s > 0 && make_stable){
			SpinFlipsolve(In, &teuk_in_HBL_integrate);
		}else{
			teuk_in_HBL_integrate(_inSolution, _inDerivative, *this, _radialPoints);
		}
	}else{
		if(_s < 0 && make_stable){
			SpinFlipsolve(Up, &teuk_up_HBL_integrate);
		}else{
			teuk_up_HBL_integrate(_upSolution, _upDerivative, *this, _radialPoints);
		}
	}
}

// Helper function for solving the Teukolsky equation. The make_stable option indicates whether or not we only 
// integrate the spin-weight solutions that are numerically stable for a given
// boundary condition
void RadialTeukolsky::TEUKsolve(BoundaryCondition bc, bool make_stable){
	if(bc == In){
		if(_s > 0 && make_stable){
			SpinFlipsolve(In, &teuk_in_TEUK_integrate);
		}else{
			teuk_in_TEUK_integrate(_inSolution, _inDerivative, *this, _radialPoints);
		}
	}else{
		if(_s < 0 && make_stable){
			SpinFlipsolve(Up, &teuk_up_TEUK_integrate);
		}else{
			teuk_up_TEUK_integrate(_upSolution, _upDerivative, *this, _radialPoints);
		}
	}
}

// Helper function for solving the Teukolsky equation using the Sasaki-Nakamura
// transformation. The make_stable option indicates whether or not we only 
// integrate the spin-weight solutions that are numerically stable for a given
// boundary condition
void RadialTeukolsky::GSNsolve(BoundaryCondition bc){
	if(_s == -2){
		if(bc == In){
			teuk_in_GSN_integrate(_inSolution, _inDerivative, *this, _radialPoints);
		}else{
			teuk_up_GSN_integrate(_upSolution, _upDerivative, *this, _radialPoints);
		}
	}else if(_s == 2){
		if(bc == In){
			SpinFlipsolve(In, &teuk_in_GSN_integrate);
		}else{
			SpinFlipsolve(Up, &teuk_up_GSN_integrate);
		}
	}
}

// Helper function that flips the sign of the spin-weight of our initial data, integrates
// the differential equation using the opposite spin-weight, then flips the spin-weight
// of the numerical solutions back to the original spin-weight using the Teukolsky-Starobinsky
// transformations
void RadialTeukolsky::SpinFlipsolve(BoundaryCondition bc, int (*func)(ComplexVector &, ComplexVector &, RadialTeukolsky&, const Vector &)){
	double lambdaCH = _lambda + _s*(_s + 1.);
	if(bc == In){
		Complex boundaryRStore = _horizonBoundarySolution.getValue();
		Complex boundaryRpStore = _horizonBoundaryDerivative.getValue();

		// flip everything to the opposite spin-weight
		Complex boundaryR, boundaryRp;
		// flip_spin_of_radial_teukolsky_TS(boundaryR, boundaryRp, In, _s, _m, _a, _omega, lambdaCH, _horizonBoundary, _horizonBoundarySolution.getValue(), _horizonBoundaryDerivative.getValue());
		// _horizonBoundarySolution = Result(boundaryR, DBL_EPSILON);
		// _horizonBoundaryDerivative = Result(boundaryRp, DBL_EPSILON);
		_lambda = flip_eigenvalue(_s, _lambda);
		_s = flip_spin(_s);
		// _horizonBoundarySolution = teuk_in_asymptotic_horizon(_a, _s, _L, _m, _omega, _lambda, _horizonBoundary);
		// _horizonBoundaryDerivative = teuk_in_derivative_asymptotic_horizon(_a, _s, _L, _m, _omega, _lambda, _horizonBoundary);
		generateRetardedBoundaryCondition(In);
		func(_inSolution, _inDerivative, *this, _radialPoints);

		// flip back everything to the original spin-weight
		_horizonBoundarySolution = Result(boundaryRStore, DBL_EPSILON);
		_horizonBoundaryDerivative = Result(boundaryRpStore, DBL_EPSILON);
		ComplexVector inSolutionTemp = _inSolution; // not sure if this copying is necessary
		ComplexVector inDerivativeTemp = _inDerivative;
		flip_spin_of_radial_teukolsky_TS(_inSolution, _inDerivative, In, _s, _m, _a, _omega, lambdaCH, _radialPoints, inSolutionTemp, inDerivativeTemp);
		_lambda = flip_eigenvalue(_s, _lambda);
		_s = flip_spin(_s);
	}else{
		Complex boundaryRStore = _infinityBoundarySolution.getValue();
		Complex boundaryRpStore = _infinityBoundaryDerivative.getValue();

		// flip everything to the opposite spin-weight
		Complex boundaryR, boundaryRp;
		// flip_spin_of_radial_teukolsky_TS(boundaryR, boundaryRp, Up, _s, _m, _a, _omega, lambdaCH, _infinityBoundary, _infinityBoundarySolution.getValue(), _infinityBoundaryDerivative.getValue());
		// _infinityBoundarySolution = Result(boundaryR, DBL_EPSILON);
		// _infinityBoundaryDerivative = Result(boundaryRp, DBL_EPSILON);

		_lambda = flip_eigenvalue(_s, _lambda);
		_s = flip_spin(_s);
		// _infinityBoundarySolution = teuk_up_asymptotic_infinity(_a, _s, _L, _m, _omega, _lambda, _infinityBoundary);
		// _infinityBoundaryDerivative = teuk_up_derivative_asymptotic_infinity(_a, _s, _L, _m, _omega, _lambda, _infinityBoundary);
		generateRetardedBoundaryCondition(Up);
		func(_upSolution, _upDerivative, *this, _radialPoints);

		// flip back everything to the original spin-weight
		_infinityBoundarySolution = Result(boundaryRStore, DBL_EPSILON);
		_infinityBoundaryDerivative = Result(boundaryRpStore, DBL_EPSILON);
		ComplexVector upSolutionTemp = _upSolution; // not sure if this copying is necessary
		ComplexVector upDerivativeTemp = _upDerivative;
		flip_spin_of_radial_teukolsky_TS(_upSolution, _upDerivative, Up, _s, _m, _a, _omega, lambdaCH, _radialPoints, upSolutionTemp, upDerivativeTemp);
		_lambda = flip_eigenvalue(_s, _lambda);
		_s = flip_spin(_s);
	}
}

// Helper function that uses asymptotic expansions to generate solutions
void RadialTeukolsky::ASYMsolve(){
	// std::cout << "(RADIALSOLVER) Asymptotic method auto selected. \n";
	teuk_ASYM_series(_inSolution, _inDerivative, _upSolution, _upDerivative, *this, _radialPoints);
}

// Helper function that uses asymptotic expansions to generate solutions for 
// a given boundary condition
void RadialTeukolsky::ASYMsolve(BoundaryCondition bc){
	// std::cout << "(RADIALSOLVER) Asymptotic method auto selected. \n";
	if(bc == In){
		teuk_in_ASYM_series(_inSolution, *this, _radialPoints);
		teuk_in_derivative_ASYM_series(_inDerivative, *this, _radialPoints);
	}else{
		teuk_up_ASYM_series(_upSolution, *this, _radialPoints);
		teuk_up_derivative_ASYM_series(_upDerivative, *this, _radialPoints);
	}
}

// Helper function that uses MST expansions to generate solutions
void RadialTeukolsky::MSTsolve(){
	// std::cout << "(RADIALSOLVER) MST method auto selected. \n";
	teuk_MST_series(_inSolution, _inDerivative, _upSolution, _upDerivative, *this, _radialPoints);
}

// Helper function that uses MST expansions to generate solutions for 
// a given boundary condition
void RadialTeukolsky::MSTsolve(BoundaryCondition bc){
	// std::cout << "(RADIALSOLVER) MST method auto selected. \n";
	if(bc == In){
		teuk_in_MST_series(_inSolution, _inDerivative, *this, _radialPoints);
	}else{
		teuk_up_MST_series(_upSolution, _upDerivative, *this, _radialPoints);
	}
}

// Helper function that gives the static solutions to the Teukolsky equation,
// which are known analytically
void RadialTeukolsky::STATICsolve(){
	// std::cout << "(RADIALSOLVER) MST method auto selected. \n";
	teuk_static(_inSolution, _inDerivative, _upSolution, _upDerivative, *this, _radialPoints);
}

// Helper function that gives the static solutions to the Teukolsky equation,
// which are known analytically, for a given boundary condition
void RadialTeukolsky::STATICsolve(BoundaryCondition bc){
	// std::cout << "(RADIALSOLVER) MST method auto selected. \n";
	teuk_static(_inSolution, _inDerivative, _upSolution, _upDerivative, *this, _radialPoints);
}

// Helper function that cycles through different solution methods and picks the
// the solution method that is the fastest and most stable depending on the 
// parameter values given to the class
void RadialTeukolsky::AUTOsolve(){
	// std::cout << "(RADIALSOLVER) AUTO method selected. \n";
	int AUTO_SUCCESS = 0;
	if(validMethodDomain(ASYM)){
		ASYMsolve();
		AUTO_SUCCESS = 1;
		if(failCheck()){ // if method failed, say auto method is not yet successful
			AUTO_SUCCESS = 0;
		}
	}

	// if(!AUTO_SUCCESS){
	// 	if(validMethodDomain(MST)){
	// 		MSTsolve();
	// 		AUTO_SUCCESS = 1;
	// 		if(failCheck()){ // if method failed, say auto method is not yet successful
	// 			AUTO_SUCCESS = 0;
	// 		}
	// 	}
	// }

	if(!AUTO_SUCCESS){
		if( _horizonBoundary == 0. || _infinityBoundary == 0. ){
			// std::cout << "(RADIALSOLVER) Numerical integration auto selected. \n";
			generateRetardedBoundaryConditions(AUTO);
		}
		HBLsolve(In, true);
		HBLsolve(Up, true);
	}
} 

// Helper function that cycles through different solution methods and picks the
// the solution method that is the fastest and most stable depending on the 
// parameter values given to the class
void RadialTeukolsky::AUTOsolve(BoundaryCondition bc){
	// std::cout << "(RADIALSOLVER) AUTO method selected. \n";
	int AUTO_SUCCESS = 0;
	if(validMethodDomain(ASYM)){
		ASYMsolve(bc);
		AUTO_SUCCESS = 1;
		if(failCheck()){ // if method failed, say auto method is not yet successful
			AUTO_SUCCESS = 0;
		}
	}

	// if(!AUTO_SUCCESS){
	// 	if(validMethodDomain(MST)){
	// 		MSTsolve(bc);
	// 		AUTO_SUCCESS = 1;
	// 		if(failCheck()){ // if method failed, say auto method is not yet successful
	// 			AUTO_SUCCESS = 0;
	// 		}
	// 	}
	// }

	if(!AUTO_SUCCESS){
		if( _horizonBoundary == 0. || _infinityBoundary == 0. ){
			// std::cout << "(RADIALSOLVER) Numerical integration auto selected. \n";
			generateRetardedBoundaryConditions(AUTO);
		}
		HBLsolve(bc, true);
	}
} 

// Generate homogeneous Teukolsky solutions using the specified method. The default method is AUTO
// which numerically integrates the hyperboloidally transformed Teukolsky equation.
// The make_stable flag is set to true by default. Therefore, for the numerical integrators will always
// integrate the Teukolsky equation using a stable choice of spin-weight and then generates the requested
// solution via the Teukolsky-Starobinsky transformations
void RadialTeukolsky::generateSolutions(SolutionMethod method, bool make_stable){
	if(std::abs(_omega) <= ZERO_FREQ_MAX){
		STATICsolve();
	}else if(method == AUTO){
		AUTOsolve();
	}else if(method == MST){
		MSTsolve();
	}else if(method == ASYM){
		ASYMsolve();
	}else{
		if( _horizonBoundary == 0. || _infinityBoundary == 0. ){
			generateRetardedBoundaryConditions(method);
		}
		if(method == HBL){
			HBLsolve(In, make_stable);
			HBLsolve(Up, make_stable);
		}else if(method == TEUK){
			TEUKsolve(In, make_stable);
			TEUKsolve(Up, make_stable);
		}else if(method == GSN && std::abs(_s) == 2){
			GSNsolve(In);
			GSNsolve(Up);
		}else{
			std::cout << "(RADIALSOLVER) ERROR: Method " << method << " is unknown or not yet implemented. Please choose another method.\n";
		}
	}
}

void RadialTeukolsky::generateSolutions(BoundaryCondition bc, SolutionMethod method, bool make_stable){
	if(std::abs(_omega) <= ZERO_FREQ_MAX){
		STATICsolve(bc);
	}else if(method == AUTO){
		AUTOsolve(bc);
	}else if(method == MST){
		MSTsolve(bc);
	}else if(method == ASYM){
		ASYMsolve(bc);
	}else{
		if( _horizonBoundary == 0. || _infinityBoundary == 0. ){
			generateRetardedBoundaryConditions(method);
		}
		if(method == HBL){
			HBLsolve(bc, make_stable);
		}else if(method == TEUK){
			TEUKsolve(bc, make_stable);
		}else if(method == GSN && std::abs(_s) == 2){
			GSNsolve(bc);
		}else{
			std::cout << "(RADIALSOLVER) ERROR: Method " << method << " is unknown or not yet implemented. Please choose another method.\n";
		}
	}
}

int RadialTeukolsky::resampleSolutions(Vector radialSamples){
	if(radialSamples.back() > _radialPoints.back() || radialSamples[0] < _radialPoints[0]){
		std::cout << "(RADIALSOLVER) ERROR: New samples must lie in the range "<<_radialPoints[0]<<" <= r <= "<<_radialPoints.back()<<" \n";
		return 0;
	}

	// if(_s == 0){
	// 	_horizonBoundary = _radialPoints[0];
	// 	_infinityBoundary = _radialPoints.back();
	//
	// 	_horizonBoundarySolution = Result(_inSolution[0], _inSolution[0]*DBL_EPSILON);
	// 	_horizonBoundaryDerivative = Result(_inDerivative[0], _inDerivative[0]*DBL_EPSILON);
	//
	// 	_infinityBoundarySolution = Result(_upSolution.back(), _upSolution.back()*DBL_EPSILON);
	// 	_infinityBoundaryDerivative = Result(_upDerivative.back(), _upDerivative.back()*DBL_EPSILON);
	// }else if(_s < 0){
	// 	_horizonBoundary = _radialPoints[0];
	// 	_infinityBoundary = _radialPoints[0];
	//
	// 	_horizonBoundarySolution = Result(_inSolution[0], _inSolution[0]*DBL_EPSILON);
	// 	_horizonBoundaryDerivative = Result(_inDerivative[0], _inDerivative[0]*DBL_EPSILON);
	//
	// 	_infinityBoundarySolution = Result(_upSolution[0], _inSolution[0]*DBL_EPSILON);
	// 	_infinityBoundaryDerivative = Result(_upDerivative[0], _upDerivative[0]*DBL_EPSILON);
	// }else{
	// 	_horizonBoundary = _radialPoints.back();
	// 	_infinityBoundary = _radialPoints.back();
	//
	// 	_horizonBoundarySolution = Result(_inSolution.back(), _inSolution.back()*DBL_EPSILON);
	// 	_horizonBoundaryDerivative = Result(_inDerivative.back(), _inDerivative.back()*DBL_EPSILON);
	//
	// 	_infinityBoundarySolution = Result(_upSolution.back(), _upSolution.back()*DBL_EPSILON);
	// 	_infinityBoundaryDerivative = Result(_upDerivative.back(), _upDerivative.back()*DBL_EPSILON);
	// }

	_horizonBoundary = _radialPoints[0];
	_infinityBoundary = _radialPoints.back();

	_horizonBoundarySolution = Result(_inSolution[0], _inSolution[0]*DBL_EPSILON);
	_horizonBoundaryDerivative = Result(_inDerivative[0], _inDerivative[0]*DBL_EPSILON);

	_infinityBoundarySolution = Result(_upSolution.back(), _upSolution.back()*DBL_EPSILON);
	_infinityBoundaryDerivative = Result(_upDerivative.back(), _upDerivative.back()*DBL_EPSILON);

	_radialPoints.resize(radialSamples.size(), 0.);
	_radialPoints = radialSamples;

	_inSolution.resize(_radialPoints.size(), 0.);
	_inDerivative.resize(_radialPoints.size(), 0.);
	teuk_in_HBL_integrate(_inSolution, _inDerivative, *this, _radialPoints);

	_upSolution.resize(_radialPoints.size(), 0.);
	_upDerivative.resize(_radialPoints.size(), 0.);
	teuk_up_HBL_integrate(_upSolution, _upDerivative, *this, _radialPoints);

	return 0;
}

// flips the spin-weight of the spin-dependent parameters, solutions, and derivatives
void RadialTeukolsky::flipSpinWeight(){
	// flip everything to the opposite spin-weight
	double lambdaCH = _lambda + _s*(_s + 1.);
	Complex boundaryR, boundaryRp;
	flip_spin_of_radial_teukolsky_TS(boundaryR, boundaryRp, In, _s, _m, _a, _omega, lambdaCH, _horizonBoundary, _horizonBoundarySolution.getValue(), _horizonBoundaryDerivative.getValue());
	_horizonBoundarySolution = Result(boundaryR, DBL_EPSILON);
	_horizonBoundaryDerivative = Result(boundaryRp, DBL_EPSILON);
	flip_spin_of_radial_teukolsky_TS(boundaryR, boundaryRp, Up, _s, _m, _a, _omega, lambdaCH, _infinityBoundary, _infinityBoundarySolution.getValue(), _infinityBoundaryDerivative.getValue());
	_infinityBoundarySolution = Result(boundaryR, DBL_EPSILON);
	_infinityBoundaryDerivative = Result(boundaryRp, DBL_EPSILON);

	ComplexVector inSolutionTemp = _inSolution;
	ComplexVector upSolutionTemp = _upSolution;
	ComplexVector inDerivativeTemp = _inDerivative;
	ComplexVector upDerivativeTemp = _upDerivative;
	flip_spin_of_radial_teukolsky_TS(_inSolution, _inDerivative, In, _s, _m, _a, _omega, lambdaCH, _radialPoints, inSolutionTemp, inDerivativeTemp);
	flip_spin_of_radial_teukolsky_TS(_upSolution, _upDerivative, Up, _s, _m, _a, _omega, lambdaCH, _radialPoints, upSolutionTemp, upDerivativeTemp);
	_lambda = flip_eigenvalue(_s, _lambda);
	_s = flip_spin(_s);
}


Vector RadialTeukolsky::getRadialPoints(){ return _radialPoints; }
double RadialTeukolsky::getRadialPoints(int pos){ return _radialPoints[pos]; }

double RadialTeukolsky::getBoundaryPoint(BoundaryCondition bc){
	if(bc == In){
		return _horizonBoundary;
	}else{
		return _infinityBoundary;
	}
}
Result RadialTeukolsky::getBoundarySolution(BoundaryCondition bc){
	if(bc == In){
		// if(_horizonBoundary == 0. || _infinityBoundary == 0. ){
		// 	generateRetardedBoundaryConditions();
		// }
		return _horizonBoundarySolution;
	}else{
		// if(_horizonBoundary == 0. || _infinityBoundary == 0. ){
		// 	generateRetardedBoundaryConditions();
		// }
		return _infinityBoundarySolution;
	}
}
Result RadialTeukolsky::getBoundaryDerivative(BoundaryCondition bc){
	if(bc == In){
		// if(_horizonBoundary == 0. || _infinityBoundary == 0. ){
		// 	generateRetardedBoundaryConditions();
		// }
		return _horizonBoundaryDerivative;
	}else{
		// if(_horizonBoundary == 0. || _infinityBoundary == 0. ){
		// 	generateRetardedBoundaryConditions();
		// }
		return _infinityBoundaryDerivative;
	}
}

ComplexVector RadialTeukolsky::getSolution(BoundaryCondition bc){
	if(bc == In){
		// if( std::abs(_inSolution[0]) == 0. && std::abs(_inDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _inSolution;
	}else{
		// if( std::abs(_upSolution[0]) == 0. && std::abs(_upDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _upSolution;
	}
}

Complex RadialTeukolsky::getSolution(BoundaryCondition bc, int pos){
	if(bc == In){
		// if( std::abs(_inSolution[0]) == 0. && std::abs(_inDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _inSolution[pos];
	}else{
		// if( std::abs(_upSolution[0]) == 0. && std::abs(_upDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _upSolution[pos];
	}
}

ComplexVector RadialTeukolsky::getDerivative(BoundaryCondition bc){
	if(bc == In){
		// if( std::abs(_inSolution[0]) == 0. && std::abs(_inDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _inDerivative;
	}else{
		// if( std::abs(_upSolution[0]) == 0. && std::abs(_upDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _upDerivative;
	}
}

Complex RadialTeukolsky::getDerivative(BoundaryCondition bc, int pos){
	if(bc == In){
		// if( std::abs(_inSolution[0]) == 0. && std::abs(_inDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _inDerivative[pos];
	}else{
		// if( std::abs(_upSolution[0]) == 0. && std::abs(_upDerivative[0]) == 0. ){
		// 	generateSolutions(bc);
		// }

		return _upDerivative[pos];
	}
}

ComplexVector RadialTeukolsky::getSecondDerivative(BoundaryCondition bc){
	ComplexVector Rpp(getDerivative(bc));
	for(size_t i = 0; i < Rpp.size(); i++){
		Rpp[i] = teuk_secondDerivative(_a, _s, _m, _omega, _lambda, _radialPoints[i], getSolution(bc)[i], getDerivative(bc)[i]);
	}
	return Rpp;
}

Complex RadialTeukolsky::getSecondDerivative(BoundaryCondition bc, int pos){
	return teuk_secondDerivative(_a, _s, _m, _omega, _lambda, _radialPoints[pos], getSolution(bc, pos), getDerivative(bc, pos));
}

Complex teuk_secondDerivative(double a, int s, int m, double omega, double lambda, double r, Complex R, Complex Rp){
	double delta = pow(r, 2) + pow(a, 2) - 2.*r;
	double K = (pow(r, 2) + pow(a, 2))*omega - m*a;
	return (-2.*(Complex(s) + 1.)*(r - 1.)*Rp - ((pow(K, 2) - 2.*I*Complex(s)*(r - 1.)*K)/delta + 4.*I*Complex(s)*omega*r - lambda)*R)/delta;
}

//*************************************************************//
//       Teukolsky-Starobinsky Identities for s = \pm 2        //
//*************************************************************//

int flip_spin(int s){
	return -s;
}
double flip_eigenvalue(int s, double lambda){
	return lambda + 2.*s;
}

// Squared norm of the Teukolsky-Starobinsky constant given in terms of (a, m, omega) 
// and Chandrasekhar's spin-invariant eigenvalue lambdaCH = lambda_s + s(s+1) for s = pm 2
double teukolsky_starobinsky_spin_2_constant(int m, double a, double omega, double lambdaCH){
	return (pow(lambdaCH, 2) + 4*m*a*omega - pow(2*a*omega, 2))*(pow(lambdaCH - 2., 2) + 36.*m*a*omega - pow(6.*a*omega, 2)) + (2.*lambdaCH - 1.)*(96.*pow(a*omega, 2) - 48.*m*a*omega) + pow(12.*omega, 2)*(1. - a*a);
}

// Squared norm of the Teukolsky-Starobinsky constant given in terms of (a, m, omega) 
// and Chandrasekhar's spin-invariant eigenvalue lambdaCH = lambda_s + s(s+1) for s = pm 1
double teukolsky_starobinsky_spin_1_constant(int m, double a, double omega, double lambdaCH){
	return lambdaCH*lambdaCH + 4*a*omega*(m - a*omega);
}

// Squared norm of the Teukolsky-Starobinsky constant given in terms of (a, m, omega) 
// and Chandrasekhar's spin-invariant eigenvalue lambdaCH = lambda_s + s(s+1)
double teukolsky_starobinsky_constant(int s, int m, double a, double omega, double lambdaCH){
	if(std::abs(s) == 2){
		return teukolsky_starobinsky_spin_2_constant(m, a, omega, lambdaCH);
	}else if(std::abs(s) == 1){
		return teukolsky_starobinsky_spin_1_constant(m, a, omega, lambdaCH);
	}else if(std::abs(s) == 0){
		return 1.;
	}else{
		std::cout << "(ERROR) s = " << s << " is not currently supported for the TS constant\n";
	}
	return 0.;
}

// overloading function to default to s=pm 2 case
double teukolsky_starobinsky_constant(int m, double a, double omega, double lambdaCH){
	return teukolsky_starobinsky_spin_2_constant(m, a, omega, lambdaCH);
}

// Real part of the Teukolsky-Starobinsky constant given in terms of (a, m, omega) 
// and Chandrasekhar's spin-invariant eigenvalue lambdaCH = lambda_s + s(s+1) for s = pm 2
double teukolsky_starobinsky_constant_D(int m, double a, double omega, double lambdaCH){
	return sqrt((pow(lambdaCH, 2) + 4*m*a*omega - pow(2*a*omega, 2))*(pow(lambdaCH - 2., 2) + 36.*m*a*omega - pow(6.*a*omega, 2))
		+ (2.*lambdaCH - 1.)*(96.*pow(a*omega, 2) - 48.*m*a*omega) - pow(12.*a*omega, 2));
}

// Complex Teukolsky-Starobinsky constant given in terms of (a, m, omega) 
// and Chandrasekhar's spin-invariant eigenvalue lambdaCH = lambda_s + s(s+1) for s = pm 2
Complex teukolsky_starobinsky_complex_constant(int j, int m, double a, double omega, double lambdaCH){
	return teukolsky_starobinsky_constant_D(m, a, omega, lambdaCH) + pow(-1., j + m)*12.*I*omega;
}

// The amplitudes of the homogeneous solutions after their spin-weights have been flipped by
// the Teukolsky-Starobinsky radial transformation for s = pm 2
Complex teukolsky_starobinsky_spin_2_amplitude(BoundaryCondition bc, int s, int m, double a, double omega, double lambdaCH){
	if(bc == In){
	if(s > 0){
		return teukolsky_starobinsky_plus_2_in(m, a, omega, lambdaCH);
	}else{
		return teukolsky_starobinsky_minus_2_in(m, a, omega, lambdaCH);
	}
	}else{
	if(s > 0){
		return teukolsky_starobinsky_plus_2_up(m, a, omega, lambdaCH);
	}else{
		return teukolsky_starobinsky_minus_2_up(m, a, omega, lambdaCH);
	}
	}
}

// The amplitudes of the homogeneous solutions after their spin-weights have been flipped by
// the Teukolsky-Starobinsky radial transformation for s = pm 1
Complex teukolsky_starobinsky_spin_1_amplitude(BoundaryCondition bc, int s, int m, double a, double omega, double lambdaCH){
  if(bc == In){
    if(s > 0){
      return teukolsky_starobinsky_plus_1_in(m, a, omega, lambdaCH);
    }else{
      return teukolsky_starobinsky_minus_1_in(m, a, omega, lambdaCH);
    }
  }else{
    if(s > 0){
      return teukolsky_starobinsky_plus_1_up(m, a, omega, lambdaCH);
    }else{
      return teukolsky_starobinsky_minus_1_up(m, a, omega, lambdaCH);
    }
  }
}

// The amplitudes of the homogeneous solutions after their spin-weights have been flipped by
// the Teukolsky-Starobinsky radial transformation for s = pm 1, 2
Complex teukolsky_starobinsky_amplitude(BoundaryCondition bc, int s, int m, double a, double omega, double lambdaCH){
	if(std::abs(s) == 0){
		return 1.;
	}else if(std::abs(s) == 1){
		return teukolsky_starobinsky_spin_1_amplitude(bc, s, m, a, omega, lambdaCH);
	}else if(std::abs(s) == 2){
		return teukolsky_starobinsky_spin_2_amplitude(bc, s, m, a, omega, lambdaCH);
	}else{
		std::cout << "(ERROR) s = " << s << " is not currently supported for the TS constant\n";
	}
	return 0.;
}

Complex teukolsky_starobinsky_plus_2_in(int m, double a, double omega, double){
  if(std::abs(omega) > 0.){
	double kappa = sqrt(1. - a*a);
	double w = (2.*omega*(1. + kappa) - m*a)/kappa;
	Complex iw = -I*w;
    return pow(2.*kappa, 4)*(iw - 1.)*(iw)*(iw + 1.)*(iw + 2.);
  }
  return 1.;
}

Complex teukolsky_starobinsky_minus_2_in(int m, double a, double omega, double lambdaCH){
  	return teukolsky_starobinsky_constant(m, a, omega, lambdaCH)/teukolsky_starobinsky_plus_2_in(m, a, omega, lambdaCH);
}

Complex teukolsky_starobinsky_plus_2_up(int m, double a, double omega, double lambdaCH){
  	return teukolsky_starobinsky_constant(m, a, omega, lambdaCH)/teukolsky_starobinsky_minus_2_up(m, a, omega, lambdaCH);
}

Complex teukolsky_starobinsky_minus_2_up(int, double, double omega, double lambdaCH){
  if(std::abs(omega) > 0.){
    return pow(2.*omega, 4);
  }
  return lambdaCH*(lambdaCH - 2.);
}

Complex teukolsky_starobinsky_plus_1_in(int m, double a, double omega, double){
  double kappa = sqrt(1. - a*a);
  double w = (2.*omega*(1. + kappa) - m*a)/kappa;
  Complex iw = -I*w;
  if(std::abs(omega) > 0.){
    return pow(2.*kappa, 2)*(iw + 1.)*(iw);
  }
  return 1.;
}

Complex teukolsky_starobinsky_minus_1_in(int m, double a, double omega, double lambdaCH){
  	return teukolsky_starobinsky_spin_1_constant(m, a, omega, lambdaCH)/teukolsky_starobinsky_plus_1_in(m, a, omega, lambdaCH);
}

Complex teukolsky_starobinsky_plus_1_up(int m, double a, double omega, double lambdaCH){
  	return teukolsky_starobinsky_spin_1_constant(m, a, omega, lambdaCH)/teukolsky_starobinsky_minus_1_up(m, a, omega, lambdaCH);
}

Complex teukolsky_starobinsky_minus_1_up(int, double, double omega, double lambdaCH){
  if(std::abs(omega) > 0.){
    return -pow(2.*omega, 2);
  }
  return lambdaCH;
}

Complex spin_minus_1_f_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=-2}
	double delta = r*r - 2.*r + a*a;
  	double Kteuk = (r*r + a*a)*omega - m*a;
	return pow(delta,-1)*(lambdaCH + 2.*I*omega*r - 2.*pow(Kteuk, 2)/delta);
}

Complex spin_minus_1_g_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=-2}
	double delta = r*r - 2.*r + a*a;
  	double Kteuk = (r*r + a*a)*omega - m*a;
	return -2.*I*Kteuk/delta;
}

Complex spin_minus_1_df_coeff_TS_static(double a, double lambdaCH, double r){ // this is lambda_{s=-2}
	double delta = r*r - 2.*r + a*a;
	return -2.*lambdaCH*(-1. + r)*pow(delta,-2);
}

Complex spin_plus_1_f_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=+2}=lambda_{s=-2}-4
	double delta = r*r - 2.*r + a*a;
	return pow(delta, 2)*std::conj(spin_minus_1_f_coeff_TS(m, a, omega, lambdaCH, r) + 2.*(r - 1.)/delta*spin_minus_1_g_coeff_TS(m, a, omega, lambdaCH, r));
}

Complex spin_plus_1_g_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=+2}=lambda_{s=-2}-4
	double delta = r*r - 2.*r + a*a;
	return pow(delta, 2)*std::conj(spin_minus_1_g_coeff_TS(m, a, omega, lambdaCH, r));
}

Complex spin_plus_1_df_coeff_TS_static(double a, double lambdaCH, double r){ // this is lambda_{s=+2}=lambda_{s=-2}-4
	double delta = r*r - 2.*r + a*a;
	return 2.*lambdaCH*(-1. + r);
}


Complex spin_minus_2_f_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=-2}
	double delta = r*r - 2.*r + a*a;
  	double Kteuk = (r*r + a*a)*omega - m*a;
	return pow(delta,-3)*(-4.*Kteuk*(-1. + r)*(Complex(0.,1.)*(2. - lambdaCH) +  8.*omega*r) - 8.*(lambdaCH - 1. + Complex(0.,3.)*omega*r)*pow(Kteuk,2)) + pow(delta,-4)*(8.*pow(Kteuk,4) + 8.*pow(Kteuk,2)*pow(-1. + r,2)) + pow(delta,-2)*((lambdaCH - 2.)*lambdaCH + 20.*Kteuk*omega + Complex(0.,4.)*omega*(3. + (lambdaCH - 2.)*r) + 12.*pow(omega,2)*pow(r,2));
}

Complex spin_minus_2_g_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=-2}
	double delta = r*r - 2.*r + a*a;
  	double Kteuk = (r*r + a*a)*omega - m*a;
	return 4.*I*((2.*pow(r - 1., 2) - lambdaCH*delta)*Kteuk + 2.*pow(Kteuk, 3) + 2.*omega*delta*(delta - r*(r - 1.)))/pow(delta, 3);
}

Complex spin_minus_2_df_coeff_TS_static(double a, double lambdaCH, double r){ // this is lambda_{s=-2}
	double delta = r*r - 2.*r + a*a;
	return -4.*lambdaCH*(lambdaCH - 2.)*(-1. + r)*pow(delta,-3);
}

Complex spin_plus_2_f_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=+2}=lambda_{s=-2}-4
	double delta = r*r - 2.*r + a*a;
	return pow(delta, 4)*std::conj(spin_minus_2_f_coeff_TS(m, a, omega, lambdaCH, r) + 4.*(r - 1.)/delta*spin_minus_2_g_coeff_TS(m, a, omega, lambdaCH, r));
}

Complex spin_plus_2_g_coeff_TS(int m, double a, double omega, double lambdaCH, double r){ // this is lambda_{s=+2}=lambda_{s=-2}-4
	double delta = r*r - 2.*r + a*a;
	return pow(delta, 4)*std::conj(spin_minus_2_g_coeff_TS(m, a, omega, lambdaCH, r));
}

Complex spin_plus_2_df_coeff_TS_static(double a, double lambdaCH, double r){ // this is lambda_{s=+2}=lambda_{s=-2}-4
	double delta = r*r - 2.*r + a*a;
	return 4.*lambdaCH*(lambdaCH - 2.)*(-1. + r)*delta;
}

Complex f_coeff_TS(int s, int m, double a, double omega, double lambdaCH, double r){
	if(std::abs(s) == 0){
		return 1.;
	}else if(std::abs(s) == 1){
		if(s < 0){
			return spin_minus_1_f_coeff_TS(m, a, omega, lambdaCH, r);
		}else{
			return spin_plus_1_f_coeff_TS(m, a, omega, lambdaCH, r);
		}
	}else if(std::abs(s) == 2){
		if(s < 0){
			return spin_minus_2_f_coeff_TS(m, a, omega, lambdaCH, r);
		}else{
			return spin_plus_2_f_coeff_TS(m, a, omega, lambdaCH, r);
		}
	}else{
		std::cout << "(ERROR) s = " << s << " is not currently supported for the TS constant\n";
	}
	return 0.;
}

Complex g_coeff_TS(int s, int m, double a, double omega, double lambdaCH, double r){
	if(std::abs(s) == 0){
		return 0.;
	}else if(std::abs(s) == 1){
		if(s < 0){
			return spin_minus_1_g_coeff_TS(m, a, omega, lambdaCH, r);
		}else{
			return spin_plus_1_g_coeff_TS(m, a, omega, lambdaCH, r);
		}
	}else if(std::abs(s) == 2){
		if(s < 0){
			return spin_minus_2_g_coeff_TS(m, a, omega, lambdaCH, r);
		}else{
			return spin_plus_2_g_coeff_TS(m, a, omega, lambdaCH, r);
		}
	}else{
		std::cout << "(ERROR) s = " << s << " is not currently supported for the TS constant\n";
	}
	return 0.;
}

Complex df_coeff_TS_static(int s, double a, double lambdaCH, double r){
	if(std::abs(s) == 0){
		return 0.;
	}else if(std::abs(s) == 1){
		if(s < 0){
			return spin_minus_1_df_coeff_TS_static(a, lambdaCH, r);
		}else{
			return spin_plus_1_df_coeff_TS_static(a, lambdaCH, r);
		}
	}else if(std::abs(s) == 2){
		if(s < 0){
			return spin_minus_2_df_coeff_TS_static(a, lambdaCH, r);
		}else{
			return spin_plus_2_df_coeff_TS_static(a, lambdaCH, r);
		}
	}else{
		std::cout << "(ERROR) s = " << s << " is not currently supported for the TS constant\n";
	}
	return 0.;
}


// void flip_spin_of_radial_teukolsky_TS(Complex &RinFlip, Complex &RupFlip, int m, double a, double omega, double lambda, double r, Complex Rin, Complex RinP, Complex Rup, Complex RupP){
//   Complex Cp2in = teukolsky_starobinsky_plus_2_in(m, a, omega, lambda);
//   Complex Cp2up = teukolsky_starobinsky_plus_2_up(m, a, omega, lambda);

//   double delta = r*r - 2.*r + a*a;
//   double Kteuk = (r*r + a*a)*omega - m*a;

//   Complex g0 = 4.*I*((2.*pow(r - 1., 2) - (lambda + 2.)*delta)*Kteuk + 2.*pow(Kteuk, 3) + 2.*omega*delta*(delta - r*(r - 1.)))/pow(delta, 3);
//   Complex f0 = pow(delta,-3)*(-4.*Kteuk*(-1. + r)*(Complex(0.,2.) - Complex(0.,1.)*(2 + lambda) +  8.*omega*r) - 8.*(1. + lambda + Complex(0.,3.)*omega*r)*pow(Kteuk,2)) + pow(delta,-4)*(8.*pow(Kteuk,4) + 8.*pow(Kteuk,2)*pow(-1. + r,2)) + pow(delta,-2)*((0. + lambda)*(2 + lambda) + 20.*Kteuk*omega + Complex(0.,4.)*omega*(3. + (0. + lambda)*r) + 12.*pow(omega,2)*pow(r,2));

//   RinFlip = (f0*Rin + g0*RinP)/Cp2in;
//   RupFlip = (f0*Rup + g0*RupP)/Cp2up;
// }

void flip_spin_of_radial_teukolsky_TS_subfunc(Complex &RFlip, Complex &RPFlip, Complex Cm, Complex Cp, int s, int m, double a, double omega, double lambdaCH, double r, Complex R, Complex RP){
	Complex g0p, f0p, g0m, f0m;
	if(std::abs(s) == 2){
		g0m = spin_minus_2_g_coeff_TS(m, a, omega, lambdaCH, r);
		f0m = spin_minus_2_f_coeff_TS(m, a, omega, lambdaCH, r);
		g0p = spin_plus_2_g_coeff_TS(m, a, omega, lambdaCH, r);
		f0p = spin_plus_2_f_coeff_TS(m, a, omega, lambdaCH, r);
	}else if(std::abs(s) == 1){
		g0m = spin_minus_1_g_coeff_TS(m, a, omega, lambdaCH, r);
		f0m = spin_minus_1_f_coeff_TS(m, a, omega, lambdaCH, r);
		g0p = spin_plus_1_g_coeff_TS(m, a, omega, lambdaCH, r);
		f0p = spin_plus_1_f_coeff_TS(m, a, omega, lambdaCH, r);
	}else{
		g0m = 0.;
		f0m = 1.;
		g0p = 0.;
		f0p = 1.;
	}

	double delta = r*r - 2.*r + a*a;
	if(s < 0){
		RFlip = (f0m*R + g0m*RP)/Cp;
		if(std::abs(g0p) > 0. && std::abs(omega) > ZERO_FREQ_MAX){
			RPFlip = (Cm*R - f0p*RFlip)/g0p;
		}else{
			Complex df0m = df_coeff_TS_static(s, a, lambdaCH, r);
			RPFlip = (df0m*R + f0m*RP)/Cp;
		}
	}else{
		RFlip = (f0p*R + g0p*RP)/Cm;
		if(std::abs(g0m) > 0. && std::abs(omega) > ZERO_FREQ_MAX){
			RPFlip = (Cp*R - f0m*RFlip)/g0m;
		}else{
			Complex df0p = df_coeff_TS_static(s, a, lambdaCH, r);
			RPFlip = (df0p*R + f0p*RP)/Cm;
		}
	}
}

void flip_spin_of_radial_teukolsky_TS(Complex &RFlip, Complex &RPFlip, BoundaryCondition bc, int s, int m, double a, double omega, double lambdaCH, double r, Complex R, Complex RP){
  	Complex Cp, Cm;
	Cp = teukolsky_starobinsky_amplitude(bc, std::abs(s), m, a, omega, lambdaCH);
	Cm = teukolsky_starobinsky_amplitude(bc, -std::abs(s), m, a, omega, lambdaCH);
	flip_spin_of_radial_teukolsky_TS_subfunc(RFlip, RPFlip, Cm, Cp, s, m, a, omega, lambdaCH, r, R, RP);
}

void flip_spin_of_radial_teukolsky_TS(ComplexVector &RFlip, ComplexVector &RPFlip, BoundaryCondition bc, int s, int m, double a, double omega, double lambdaCH, Vector r, ComplexVector R, ComplexVector RP){
  	Complex Cp, Cm;
	Cp = teukolsky_starobinsky_amplitude(bc, std::abs(s), m, a, omega, lambdaCH);
	Cm = teukolsky_starobinsky_amplitude(bc, -std::abs(s), m, a, omega, lambdaCH);
	for(size_t i = 0; i < r.size(); i++){
		flip_spin_of_radial_teukolsky_TS_subfunc(RFlip[i], RPFlip[i], Cm, Cp, s, m, a, omega, lambdaCH, r[i], R[i], RP[i]);
	}
}

// void flip_spin_of_radial_teukolsky_TS(Complex &RinFlip, Complex &RinPFlip, Complex &RupFlip, Complex &RupPFlip, int m, double a, double omega, double lambda, double r, Complex Rin, Complex RinP, Complex Rup, Complex RupP){
//   Complex Cp2in = teukolsky_starobinsky_amplitude(In, 2, m, a, omega, lambda);
//   Complex Cp2up = teukolsky_starobinsky_amplitude(Up, 2, m, a, omega, lambda);
//   Complex Cm2in = teukolsky_starobinsky_amplitude(In, -2, m, a, omega, lambda);
//   Complex Cm2up = teukolsky_starobinsky_amplitude(Up, -2, m, a, omega, lambda);

//   double delta = r*r - 2.*r + a*a;
//   double Kteuk = (r*r + a*a)*omega - m*a;

//   Complex g0 = 4.*I*((2.*pow(r - 1., 2) - (lambda + 2.)*delta)*Kteuk + 2.*pow(Kteuk, 3) + 2.*omega*delta*(delta - r*(r - 1.)))/pow(delta, 3);
//   Complex f0 = pow(delta,-3)*(-4.*Kteuk*(-1. + r)*(Complex(0.,2.) - Complex(0.,1.)*(2 + lambda) +  8.*omega*r) - 8.*(1. + lambda + Complex(0.,3.)*omega*r)*pow(Kteuk,2)) + pow(delta,-4)*(8.*pow(Kteuk,4) + 8.*pow(Kteuk,2)*pow(-1. + r,2)) + pow(delta,-2)*((0. + lambda)*(2 + lambda) + 20.*Kteuk*omega + Complex(0.,4.)*omega*(3. + (0. + lambda)*r) + 12.*pow(omega,2)*pow(r,2));

//   RinFlip = (f0*Rin + g0*RinP)/Cp2in;
//   RupFlip = (f0*Rup + g0*RupP)/Cp2up;

//   if(std::abs(g0) > 0.){
//     RinPFlip = (Cm2in*pow(delta, -4)*Rin - std::conj(f0 + 4.*(r - 1.)/delta*g0)*RinFlip)/std::conj(g0);
//     RupPFlip = (Cm2up*pow(delta, -4)*Rup - std::conj(f0 + 4.*(r - 1.)/delta*g0)*RupFlip)/std::conj(g0);
//   }else{
//     Complex df0 = -4.*lambda*(2. + lambda)*(-1. + r)*pow(delta,-3);
//     RinPFlip = (df0*Rin + f0*RinP)/Cp2in;
//     RupPFlip = (df0*Rup + f0*RupP)/Cp2up;
//   }
// }

// void flip_spin_of_radial_teukolsky(ComplexVector &RinFlip, ComplexVector &RinPFlip, ComplexVector &RupFlip, ComplexVector &RupPFlip, int m, double a, double omega, double lambda, Vector rVec, ComplexVector RinVec, ComplexVector RinPVec, ComplexVector RupVec, ComplexVector RupPVec){

//   double r = rVec[0];
//   double delta = r*r - 2.*r + a*a;
//   Complex Rin = RinVec[0];
//   Complex dRin = RinPVec[0];
//   Complex Rup = RupVec[0];
//   Complex dRup = RupPVec[0];
//   Complex RinFlip0, dRinFlip0, RupFlip0, dRupFlip0;

//   flip_spin_of_radial_teukolsky_TS(RinFlip0, dRinFlip0, RupFlip0, dRupFlip0, m, a, omega, lambda, r, Rin, dRin, Rup, dRup);

//   RinFlip[0] = RinFlip0;
//   RupFlip[0] = RupFlip0;
//   RinPFlip[0] = dRinFlip0;
//   RupPFlip[0] = dRupFlip0;

//   Complex RinConj = std::conj(Rin);
//   Complex RupConj = std::conj(Rup);
//   Complex dRinConj = std::conj(dRin);
//   Complex dRupConj = std::conj(dRup);

//   Complex DeltaRinFlip = pow(delta, 2)*RinFlip0;
//   Complex dDeltaRinFlip = pow(delta, 2)*dRinFlip0 + 2*(r - 1.)*delta*RinFlip0;
//   Complex DeltaRupFlip = pow(delta, 2)*RupFlip0;
//   Complex dDeltaRupFlip = pow(delta, 2)*dRupFlip0 + 2*(r - 1.)*delta*RupFlip0;

//   Complex wronskian = (RinConj*dRupConj - dRinConj*RupConj);

//   Complex Ain = (DeltaRinFlip*dRupConj - dDeltaRinFlip*RupConj)/wronskian;
//   Complex Aup = (DeltaRupFlip*dRupConj - dDeltaRupFlip*RupConj)/wronskian;
//   Complex Bin = -(DeltaRinFlip*dRinConj - dDeltaRinFlip*RinConj)/wronskian;
//   Complex Bup = -(DeltaRupFlip*dRinConj - dDeltaRupFlip*RinConj)/wronskian;

//   // std::cout << "CmInc = " << std::conj(1./Aup) << "\n";
//   // std::cout << "CmRef = " << -std::conj(Bup/Aup) << "\n";
//   // std::cout << "CpInc = " << std::conj(1./Bin) << "\n";
//   // std::cout << "CpRef = " << -std::conj(Ain/Bin) << "\n";

//   for(size_t i = 1; i < rVec.size(); i++){
//     r = rVec[i];
//     delta = r*r - 2.*r + a*a;

//     Rin = RinVec[i];
//     dRin = RinPVec[i];
//     Rup = RupVec[i];
//     dRup = RupVec[i];

//     RinConj = std::conj(Rin);
//     dRinConj = std::conj(dRin);
//     RupConj = std::conj(Rup);
//     dRupConj = std::conj(dRup);

//     RinFlip[i] = pow(delta, -2)*(Ain*RinConj + Bin*RupConj);
//     RupFlip[i] = pow(delta, -2)*(Aup*RinConj + Bup*RupConj);

//     // flip_spin_of_radial_teukolsky_TS(RinFlip1, RupFlip1, m, a, omega, lambda, r, Rin, RinP, Rup, RupP);
//     // std::cout << "r["<<i<<"] = " << r << ", RinFlip error = " << 1. - RinFlip[i]/RinFlip1 << "\n";
//     // std::cout << "r["<<i<<"] = " << r << ", RupFlip error = " << 1. - RupFlip[i]/RupFlip1 << "\n";

//     RinPFlip[i] = pow(delta, -2)*(Ain*dRinConj + Bin*dRupConj) - 4.*(r - 1.)*pow(delta, -3)*(Ain*RinConj + Bin*RupConj);
//     RupPFlip[i] = pow(delta, -2)*(Aup*dRinConj + Bup*dRupConj) - 4.*(r - 1.)*pow(delta, -3)*(Aup*RinConj + Bup*RupConj);
//   }
// }


//*************************************************************//
//*************************************************************//
// Different methods for solving the radial Teukolsky equation //
//*************************************************************//
//*************************************************************//

//*************************************************************//
// 						MST	expansions						   //
//*************************************************************//

// Use the semi-analytic expansions of Mano, Suzuki, and Takasugi to calculate radial Teukolsky solutions
int teuk_MST_series_boundary(Result &Rin, Result &RinP, Result &Rup, Result &RupP, RadialTeukolsky &teuk, const double &rIn, const double &rUp){
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex inTransmissionAmplitude = mstWorkspace.getAmplitude(In, Transmission);
	Complex upTransmissionAmplitude = mstWorkspace.getAmplitude(Up, Transmission);
	Rin = mstWorkspace.getSolution(In, rIn)/inTransmissionAmplitude;
	RinP = mstWorkspace.getDerivative(In, rIn)/inTransmissionAmplitude;
	Rup = mstWorkspace.getSolution(Up, rUp)/upTransmissionAmplitude;
	RupP = mstWorkspace.getDerivative(Up, rUp)/upTransmissionAmplitude;

	return 0;
}

int teuk_in_MST_series_boundary(Result &Rin, Result &RinP, RadialTeukolsky &teuk, const double &rIn){
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex inTransmissionAmplitude = mstWorkspace.getAmplitude(In, Transmission);
	Rin = mstWorkspace.getSolution(In, rIn)/inTransmissionAmplitude;
	RinP = mstWorkspace.getDerivative(In, rIn)/inTransmissionAmplitude;

	return 0;
}

int teuk_up_MST_series_boundary(Result &Rup, Result &RupP, RadialTeukolsky &teuk, const double &rUp){
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex upTransmissionAmplitude = mstWorkspace.getAmplitude(Up, Transmission);
	Rup = mstWorkspace.getSolution(Up, rUp)/upTransmissionAmplitude;
	RupP = mstWorkspace.getDerivative(Up, rUp)/upTransmissionAmplitude;

	return 0;
}

int teuk_MST_series(ComplexVector &Rin, ComplexVector &RinP, ComplexVector &Rup, ComplexVector &RupP, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = Rin.size();
	nmax = r.size() < nmax ? r.size() : nmax;
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex inTransmissionAmplitude = mstWorkspace.getAmplitude(In, Transmission);
	Complex upTransmissionAmplitude = mstWorkspace.getAmplitude(Up, Transmission);
	for(size_t n = 0; n < nmax; n++){
		Rin[n] = mstWorkspace.getSolution(In, r[n]).getValue()/inTransmissionAmplitude;
		RinP[n] = mstWorkspace.getDerivative(In, r[n]).getValue()/inTransmissionAmplitude;
		Rup[n] = mstWorkspace.getSolution(Up, r[n]).getValue()/upTransmissionAmplitude;
		RupP[n] = mstWorkspace.getDerivative(Up, r[n]).getValue()/upTransmissionAmplitude;
	}

	return 0;
}

int teuk_in_MST_series(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex transmissionAmplitude = mstWorkspace.getAmplitude(In, Transmission);
	for(size_t n = 0; n < nmax; n++){
		R[n] = mstWorkspace.getSolution(In, r[n]).getValue()/transmissionAmplitude;
		Rp[n] = mstWorkspace.getDerivative(In, r[n]).getValue()/transmissionAmplitude;
	}

	return 0;
}

int teuk_in_MST_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex transmissionAmplitude = mstWorkspace.getAmplitude(In, Transmission);
	for(size_t n = 0; n < nmax; n++){
		R[n] = mstWorkspace.getSolution(In, r[n]).getValue()/transmissionAmplitude;
	}

	return 0;
}
int teuk_in_derivative_MST_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = Rp.size();
	nmax = r.size() < nmax ? r.size() : nmax;
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex transmissionAmplitude = mstWorkspace.getAmplitude(In, Transmission);
	for(size_t n = 0; n < nmax; n++){
		Rp[n] = mstWorkspace.getDerivative(In, r[n]).getValue()/transmissionAmplitude;
	}

	return 0;
}
int teuk_up_MST_series(ComplexVector &R, ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex transmissionAmplitude = mstWorkspace.getAmplitude(Up, Transmission);
	for(size_t n = 0; n < nmax; n++){
		R[n] = mstWorkspace.getSolution(Up, r[n]).getValue()/transmissionAmplitude;
		Rp[n] = mstWorkspace.getDerivative(Up, r[n]).getValue()/transmissionAmplitude;
	}

	return 0;
}
int teuk_up_MST_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex transmissionAmplitude = mstWorkspace.getAmplitude(Up, Transmission);
	for(size_t n = 0; n < nmax; n++){
		R[n] = mstWorkspace.getSolution(Up, r[n]).getValue()/transmissionAmplitude;
	}

	return 0;
}
int teuk_up_derivative_MST_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = Rp.size();
	nmax = r.size() < nmax ? r.size() : nmax;
	MstSeriesWorkspace mstWorkspace(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), 2.*teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue());

	Complex transmissionAmplitude = mstWorkspace.getAmplitude(Up, Transmission);
	for(size_t n = 0; n < nmax; n++){
		Rp[n] = mstWorkspace.getDerivative(Up, r[n]).getValue()/transmissionAmplitude;
	}

	return 0;
}

//*************************************************************//
// 				Expansions near singular points				   //
//*************************************************************//

int teuk_ASYM_series(ComplexVector &Rin, ComplexVector &RinP, ComplexVector &Rup, ComplexVector &RupP, RadialTeukolsky &teuk, Vector &r){
	int check = teuk_in_ASYM_series(Rin, teuk, r);
	check += teuk_in_derivative_ASYM_series(RinP, teuk, r);
	check += teuk_up_ASYM_series(Rup, teuk, r);
	check += teuk_up_derivative_ASYM_series(RupP, teuk, r);

	return check;
}

Result teuk_in_ASYM_series(RadialTeukolsky &teuk, const double &r){
	return teuk_in_asymptotic_horizon(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue(), r);
}

int teuk_in_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_in_ASYM_series(teuk, r);
	return 0;
}

int teuk_in_ASYM_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		R[n] = teuk_in_ASYM_series(teuk, r[n]).getValue();
	}

	return 0;
}

Result teuk_in_derivative_ASYM_series(RadialTeukolsky &teuk, const double &r){
	return teuk_in_derivative_asymptotic_horizon(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue(), r);
}

int teuk_in_derivative_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_in_derivative_ASYM_series(teuk, r);
	return 0;
}

int teuk_in_derivative_ASYM_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = Rp.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		Rp[n] = teuk_in_derivative_ASYM_series(teuk, r[n]).getValue();
	}

	return 0;
}

Result teuk_up_ASYM_series(RadialTeukolsky &teuk, const double &r){
	return teuk_up_asymptotic_infinity(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue(), r);
}

int teuk_up_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_up_ASYM_series(teuk, r);
	return 0;
}

int teuk_up_ASYM_series(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		R[n] = teuk_up_ASYM_series(teuk, r[n]).getValue();
	}

	return 0;
}

Result teuk_up_derivative_ASYM_series(RadialTeukolsky &teuk, const double &r){
	return teuk_up_derivative_asymptotic_infinity(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), teuk.getModeFrequency(), teuk.getSpinWeightedSpheroidalEigenvalue(), r);
}

int teuk_up_derivative_ASYM_series(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_up_derivative_ASYM_series(teuk, r);
	return 0;
}

int teuk_up_derivative_ASYM_series(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = Rp.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		Rp[n] = teuk_up_derivative_ASYM_series(teuk, r[n]).getValue();
	}

	return 0;
}

//*************************************************************//
// 				Static solutions			   //
//*************************************************************//

int teuk_static(ComplexVector &Rin, ComplexVector &RinP, ComplexVector &Rup, ComplexVector &RupP, RadialTeukolsky &teuk, Vector &r){
	int check = teuk_in_static(Rin, teuk, r);
	check += teuk_in_derivative_static(RinP, teuk, r);
	check += teuk_up_static(Rup, teuk, r);
	check += teuk_up_derivative_static(RupP, teuk, r);

	return check;
}

Result teuk_in_static(RadialTeukolsky &teuk, const double &r){
	return teuk_in_static_solution(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), r);
}

int teuk_in_static(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_in_static(teuk, r);
	return 0;
}

int teuk_in_static(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		R[n] = teuk_in_static(teuk, r[n]).getValue();
	}

	return 0;
}

Result teuk_in_derivative_static(RadialTeukolsky &teuk, const double &r){
	return teuk_in_static_derivative(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), r);
}

int teuk_in_derivative_static(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_in_derivative_static(teuk, r);
	return 0;
}

int teuk_in_derivative_static(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = Rp.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		Rp[n] = teuk_in_derivative_static(teuk, r[n]).getValue();
	}

	return 0;
}

Result teuk_up_static(RadialTeukolsky &teuk, const double &r){
	return teuk_up_static_solution(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), r);
}

int teuk_up_static(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_up_static(teuk, r);
	return 0;
}

int teuk_up_static(ComplexVector &R, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = R.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		R[n] = teuk_up_static(teuk, r[n]).getValue();
	}

	return 0;
}

Result teuk_up_derivative_static(RadialTeukolsky &teuk, const double &r){
	return teuk_up_static_derivative(teuk.getBlackHoleSpin(), teuk.getSpinWeight(), teuk.getSpheroidalModeNumber(), teuk.getAzimuthalModeNumber(), r);
}

int teuk_up_derivative_static(Result &R, RadialTeukolsky &teuk, const double &r){
	R = teuk_up_derivative_static(teuk, r);
	return 0;
}

int teuk_up_derivative_static(ComplexVector &Rp, RadialTeukolsky &teuk, const Vector &r){
	size_t nmax = Rp.size();
	nmax = r.size() < nmax ? r.size() : nmax;

	for(size_t n = 0; n < nmax; n++){
		Rp[n] = teuk_up_derivative_static(teuk, r[n]).getValue();
	}

	return 0;
}

Result teuk_in_static_solution(const double &a, const int &s, const int &l, const int &m, const double &r){
	double kappa = sqrt(1. - a*a);
	double tau = -(m*a)/kappa;
	double x = (1. + kappa - r)/(2.*kappa);
	if(std::abs(tau) > 0 || s <= 0){
		return pow(-2.*kappa*x, -double(s) - 0.5*I*tau)*pow(2.*kappa*(1. - x), -0.5*I*tau)*hypergeo_2F1(double(-l) - I*tau, double(l+1) - I*tau, double(1 - s) - I*tau, x)/cgamma(double(1 - s) - I*tau);
	}else{
		return pow(2.*kappa, -s)*phammer(double(l - s + 1), 2*s)/double(factorial(s))*hypergeo_2F1(double(-l + s), double(l+1+s), double(1 + s), x);
	}
}

Result teuk_in_static_derivative(const double &a, const int &s, const int &l, const int &m, const double &r){
	double kappa = sqrt(1. - a*a);
	double x = (1. + kappa - r)/(2.*kappa);
	double tau = -(m*a)/kappa;
	if(std::abs(tau) > 0 || s <= 0){
		Result f1 = hypergeo_2F1(double(-l) - I*tau, double(l+1) - I*tau, double(1 - s) - I*tau, x);
		Result df1 = dhypergeo_2F1(double(-l) - I*tau, double(l+1) - I*tau, double(1 - s) - I*tau, x);
		return pow(-2.*kappa*x, -double(s) - 0.5*I*tau)*pow(2.*kappa*(1. - x), -0.5*I*tau)*((2.*double(s)*(1. - x) - (2.*x - 1.)*I*tau)*f1 - 2.*(1. - x)*x*df1)/(cgamma(double(1 - s) - I*tau)*4.*kappa*x*(1. - x));
	}else{
		return -pow(2.*kappa, -s - 1)*phammer(double(l - s + 1), 2*s)/double(factorial(s))*dhypergeo_2F1(double(-l + s), double(l+1+s), double(1 + s), x);
	}
}

Result teuk_up_static_solution(const double &a, const int &s, const int &l, const int &m, const double &r){
	double kappa = sqrt(1. - a*a);
	double x = (1. + kappa - r)/(2.*kappa);
	double tau = -(m*a)/kappa;
	return pow(-2.*kappa*x, -double(s) - 0.5*I*tau)*pow(2.*kappa*(1. - x), -double(l + 1) + 0.5*I*tau)*hypergeo_2F1(double(l + 1) - I*tau, double(l + 1 - s), double(2*l + 2), 1./(1. - x));
}

Result teuk_up_static_derivative(const double &a, const int &s, const int &l, const int &m, const double &r){
	double kappa = sqrt(1. - a*a);
	double x = (1. + kappa - r)/(2.*kappa);
	double tau = -(m*a)/kappa;
	return pow(-2.*kappa*x, -double(s) - 0.5*I*tau)*pow(2.*kappa*(1. - x), -double(l + 1) + 0.5*I*tau)*((x - 1.)*(-2.*double(s) + 2.*double(l + s + 1)*x - I*tau)*hypergeo_2F1(double(l + 1) - I*tau, double(l + 1 - s), double(2*l + 2), 1./(1. - x)) - 2.*x*dhypergeo_2F1(double(l + 1) - I*tau, double(l + 1 - s), double(2*l + 2), 1./(1. - x)))/(4.*kappa*x*pow(1. - x, 2));
}

//*************************************************************//
// 						Teukolsky ODE solvers 				   //
//*************************************************************//


//typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_vector > error_stepper_type;
//typedef boost::numeric::odeint::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
struct streaming_observer{
     std::ostream &m_out;
     streaming_observer( std::ostream &out ) : m_out( out ) {}

     void operator()( const state_type &x , double t ) const
     {
          m_out << t;
          for( size_t i=0 ; i < x.size() ; ++i )
              m_out << "\t" << x[i];
          m_out << "\n";
     }
};

template <typename ODE_FUNC>
int teuk_integrate_boost(ComplexVector &psi, ComplexVector &dpsidr, ODE_FUNC sys, state_type psi0, const double r0, const Vector &r){
	double r_step = r0;
	double first_step = r[0];
	double initial_step_size = TEUK_ODE_STEP_SIZE_INIT;
	if( first_step < r_step ){
		initial_step_size *= -1.;
	}

	int stepNum = r.size();

	double abs_err = 0., rel_err = TEUK_ODE_REL_ERR;

	//Explicit methods

	// Runga-kutta with error control
	typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_type> error_stepper_type;
	auto c_stepper = boost::numeric::odeint::make_controlled<error_stepper_type>(abs_err, rel_err);

	// Bulirsch-stoer with error control
	//boost::numeric::odeint::bulirsch_stoer_dense_out< state_type > c_stepper( abs_err , rel_err , a_x , a_dxdt );

	// Integrate equations
	state_type psi_step = psi0;
	int steps = boost::numeric::odeint::integrate_adaptive(c_stepper, sys, psi_step, r_step, r[0], initial_step_size);
	printf("Took %d steps from r_step = %f, rp = %f\n", steps, r_step, r[0]);
	psi[0] = psi_step[0] + I*psi_step[1];
	dpsidr[0] = psi_step[2] + I*psi_step[3];
	for( int i = 1; i < stepNum; i++ ){
	    steps = boost::numeric::odeint::integrate_adaptive(c_stepper, sys, psi_step, r[i-1], r[i], initial_step_size);
		psi[i] = psi_step[0] + I*psi_step[1];
		dpsidr[i] = psi_step[2] + I*psi_step[3];
	}

	return 0;
}

template <typename ODE_FUNC>
int teuk_integrate_boost_implicit(ComplexVector &psi, ComplexVector &dpsidr, ODE_FUNC sys, ODE_FUNC jac, state_type psi0, const double r0, const Vector &r){
	double r_step = r0;
	double first_step = r[0];
	double initial_step_size = TEUK_ODE_STEP_SIZE_INIT;
	if( first_step < r_step ){
		initial_step_size *= -1.;
	}
	int stepNum = r.size();

	double abs_err = TEUK_ODE_REL_ERR, rel_err = TEUK_ODE_REL_ERR;
	// Implicit methods

	//rosenbrock4 with error control
	typedef boost::numeric::odeint::rosenbrock4<double> error_stepper_type;
	auto c_stepper = boost::numeric::odeint::make_controlled<error_stepper_type>(abs_err, rel_err);

	// Integrate equations
	vector_type psi_step(4);
	for(int i = 0; i < 4; i++){
		psi_step[i] = psi0[i];
	}
	int steps = boost::numeric::odeint::integrate_adaptive(c_stepper, std::make_pair(sys, jac), psi_step, r_step, r[0], initial_step_size);

	psi[0] = psi_step[0] + I*psi_step[1];
	dpsidr[0] = psi_step[2] + I*psi_step[3];
	for( int i = 1; i < stepNum; i++ ){
	    steps = boost::numeric::odeint::integrate_adaptive(c_stepper, std::make_pair(sys, jac), psi_step, r_step, r[0], initial_step_size);
		psi[i] = psi_step[0] + I*psi_step[1];
		dpsidr[i] = psi_step[2] + I*psi_step[3];
	}

	return 0;
}

int teuk_jac_null_gsl(double, const double*, double *, double*, void*){
	return GSL_SUCCESS;
}

int teuk_integrate_gsl(ComplexVector &psi, ComplexVector &dpsidr, int (*sys)(double, const double*, double*, void*), state_type psi0, const double r0, const Vector &r, void *params){
	size_t dim = 4;
	double Psi[dim];
	for(size_t i = 0; i < dim; i++){
		Psi[i] = psi0[i];
	}

	double h = TEUK_ODE_STEP_SIZE_INIT*sqrt(Psi[0]*Psi[0] + Psi[1]*Psi[1])/sqrt(Psi[2]*Psi[2] + Psi[3]*Psi[3]);
	// double h = TEUK_ODE_STEP_SIZE_INIT;
	double rTest = r[0];
	if(std::abs(rTest - r0) <= 1.e-14){
		rTest = r[1];
	}
	if(h > std::abs(r0 - rTest)){
		h = (1.e-6)*std::abs(r0 - rTest);
	}
	if(rTest < r0){
		h *= -1;
	}
	// std::cout << "(RADIALSOLVER) Initial step size = " << h << "\n";

	gsl_odeiv2_system gsl_sys = {sys, NULL, dim, params};

	double abs_error = sqrt(Psi[0]*Psi[0] + Psi[1]*Psi[1])*TEUK_ODE_REL_ERR*(1.e-5);
	const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
	gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(abs_error, TEUK_ODE_REL_ERR);
	gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);

	int stepNum = r.size(), status;
	int stepCount = 0;
	double ri = r0;
	for(int i = 0; i < stepNum; i++){
		if( h > 0 ){
			while(ri < r[i]){
				status = gsl_odeiv2_evolve_apply(e, c, s, &gsl_sys, &ri, r[i], &h, Psi);
				// std::cout << "(RADIALSOLVER) Step size = " << h << ". Psi = "<<Psi[0] + I*Psi[1]<<" \n";
				stepCount += 1;
			}
		}else{
			while(ri > r[i]){
				status = gsl_odeiv2_evolve_apply(e, c, s, &gsl_sys, &ri, r[i], &h, Psi);
				// std::cout << "(RADIALSOLVER) Step size = " << h << ". Psi = "<<Psi[0] + I*Psi[1]<<" \n";
				stepCount += 1;
			}
		}
		psi[i] = Psi[0] + I*Psi[1];
		dpsidr[i] = Psi[2] + I*Psi[3];
	}
	// std::cout << "Step count = " << stepCount << "\n";

	gsl_odeiv2_evolve_free(e);
	gsl_odeiv2_control_free(c);
	gsl_odeiv2_step_free(s);

	return GSL_SUCCESS;
}

int teuk_integrate_gsl(ComplexVector &psi, ComplexVector &dpsidr, int (*sys)(double, const double*, double*, void*), int (*jac)(double, const double*, double*, double*, void*), state_type psi0, const double r0, const Vector &r, void *params){
	size_t dim = 4;
	double Psi[dim];
	for(size_t i = 0; i < dim; i++){
		Psi[i] = psi0[i];
	}

	double h = TEUK_ODE_STEP_SIZE_INIT/sqrt(Psi[2]*Psi[2] + Psi[3]*Psi[3]);
	double rTest = r[0];
	if(std::abs(rTest - r0) <= 1.e-14){
		rTest = r[1];
	}
	if(h > std::abs(r0 - rTest)){
		h = 0.1*std::abs(r0 - rTest);
	}
	if(rTest < r0){
		h *= -1;
	}
	//std::cout << "(RADIALSOLVER) Initial step size = " << initialStep << "\n";

	gsl_odeiv2_system gsl_sys = {sys, jac, dim, params};
	double abs_error = sqrt(Psi[0]*Psi[0] + Psi[1]*Psi[1])*TEUK_ODE_REL_ERR*(1.e-5);
	abs_error = 0.;

	const gsl_odeiv2_step_type *T = gsl_odeiv2_step_bsimp;
	gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
	gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(abs_error, TEUK_ODE_REL_ERR);
	gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);


	int stepNum = r.size();
	double ri = r0;
	int status = GSL_SUCCESS;
	for(int i = 0; i < stepNum; i++){
		if( h > 0 ){
			while(ri < r[i] && status == GSL_SUCCESS){
				status = gsl_odeiv2_evolve_apply(e, c, s, &gsl_sys, &ri, r[i], &h, Psi);
			}
		}else{
			while(ri > r[i] && status == GSL_SUCCESS){
				status = gsl_odeiv2_evolve_apply(e, c, s, &gsl_sys, &ri, r[i], &h, Psi);
			}
		}
		if(status != GSL_SUCCESS){
			printf ("GSL error: %s\n", gsl_strerror(status));
		}
		psi[i] = Psi[0] + I*Psi[1];
		dpsidr[i] = Psi[2] + I*Psi[3];
	}

	gsl_odeiv2_evolve_free (e);
	gsl_odeiv2_control_free (c);
	gsl_odeiv2_step_free (s);

	return GSL_SUCCESS;
};


//*************************************************************//
// 				Teukolsky equation in BoyerLindquist 		   //
//*************************************************************//

void teuk_spin_minus_transform_solution(Complex &R, double r, hbl_parameters params){
	double delta = pow(r, 2) - 2.*r + pow(params.a, 2);
	R = std::conj(R)*pow(delta, params.s);
}

void teuk_spin_minus_transform_derivative(Complex &Rp, Complex R, double r, hbl_parameters params){
	double delta = pow(r, 2) - 2.*r + pow(params.a, 2);
	Rp = Complex(2*params.s)*(r - 1.)*std::conj(R)*pow(delta, params.s - 1) + std::conj(Rp)*pow(delta, params.s);
}

void teuk_spin_minus_transform_parameters(hbl_parameters &params){
	params.la += 2.*params.s;
	params.s *= -1;
	params.H *= -1;
}

int teuk_in_TEUK_integrate(ComplexVector &R_return, ComplexVector &Rp_return, RadialTeukolsky &teuk, const Vector &r){
	hbl_parameters params = {
		.a = teuk.getBlackHoleSpin(),
		.s = teuk.getSpinWeight(),
		.m = teuk.getAzimuthalModeNumber(),
		.om = teuk.getModeFrequency(),
		.la = teuk.getSpinWeightedSpheroidalEigenvalue(),
		.H = 0
	};

	Complex R = teuk.getBoundarySolution(In).getValue();
	Complex Rp = teuk.getBoundaryDerivative(In).getValue();
	double r0 = teuk.getBoundaryPoint(In);

	double norm = std::abs(R);
	R = R/norm;
	Rp = Rp/norm;

	int stepNum = r.size();
	Vector rReverse = r;
	if(r0 >= r.back()){
		std::reverse(rReverse.begin(),rReverse.end());
	}

	// std::cout << std::setprecision(15);
// 	std::cout << "(RADIALSOLVER) Rin(r = "<< r0 << ") = " << R << "\n";
// 	std::cout << "(RADIALSOLVER) Rin'(r = "<< r0 << ") = " << Rp << "\n";

	state_type R0 {std::real(R), std::imag(R), std::real(Rp), std::imag(Rp)};
	if(BOOST_INTEGRATE){
		teuk_blc sys(params);
		teuk_integrate_boost(R_return, Rp_return, sys, R0, r0, rReverse);
	}else{
		teuk_integrate_gsl(R_return, Rp_return, &teuk_blc_gsl, R0, r0, rReverse, &params);
	}

	if(r0 >= r.back()){
		std::reverse(R_return.begin(),R_return.end());
		std::reverse(Rp_return.begin(),Rp_return.end());
	}

	for(int i = 0; i < stepNum; i++){
		Rp_return[i] = norm*Rp_return[i];
		R_return[i] = norm*R_return[i];
	}

	// std::cout << "(RADIALSOLVER) Rin(r = "<< r[r.size()-1] << ") = " << R_return[r.size()-1] << "\n";
// 	std::cout << "(RADIALSOLVER) Rin'(r = "<< r[r.size()-1] << ") = " << Rp_return[r.size()-1] << "\n";

	return 1;
}

int teuk_up_TEUK_integrate(ComplexVector &R_return, ComplexVector &Rp_return, RadialTeukolsky &teuk, const Vector &r){
	hbl_parameters params = {
		.a = teuk.getBlackHoleSpin(),
		.s = teuk.getSpinWeight(),
		.m = teuk.getAzimuthalModeNumber(),
		.om = teuk.getModeFrequency(),
		.la = teuk.getSpinWeightedSpheroidalEigenvalue(),
		.H = 0
	};

	Complex R = teuk.getBoundarySolution(Up).getValue();
	Complex Rp = teuk.getBoundaryDerivative(Up).getValue();
	double r0 = teuk.getBoundaryPoint(Up);

	// std::cout << std::setprecision(15);
	// std::cout << "(RADIALSOLVER) Rup(r = "<< r0 << ") = " << R << "\n";
	// std::cout << "(RADIALSOLVER) Rup'(r = "<< r0 << ") = " << Rp << "\n";

	// if(teuk.getSpinWeight() < 0){
	// 	teuk_spin_minus_transform_derivative(Rp, R, r0, params);
	// 	teuk_spin_minus_transform_solution(R, r0, params);
	// 	teuk_spin_minus_transform_parameters(params);
	// }
	double norm = std::abs(R);
	R = R/norm;
	Rp = Rp/norm;

	int stepNum = r.size();
	Vector rReverse = r;
	if(r0 >= r.back()){
		std::reverse(rReverse.begin(),rReverse.end());
	}

	state_type R0 {std::real(R), std::imag(R), std::real(Rp), std::imag(Rp)};
	if(BOOST_INTEGRATE){
		teuk_blc sys(params);
		teuk_integrate_boost(R_return, Rp_return, sys, R0, r0, rReverse);
	}else{
		teuk_integrate_gsl(R_return, Rp_return, &teuk_blc_gsl, R0, r0, rReverse, &params);
	}
	if(r0 >= r.back()){
		std::reverse(R_return.begin(),R_return.end());
		std::reverse(Rp_return.begin(),Rp_return.end());
	}

	for(int i = 0; i < stepNum; i++){
		Rp_return[i] = norm*Rp_return[i];
		R_return[i] = norm*R_return[i];
	}

	// std::cout << "(RADIALSOLVER) Rup(r = "<< r[0] << ") = " << R_return[0] << "\n";
	// std::cout << "(RADIALSOLVER) Rup'(r = "<< r[0] << ") = " << Rp_return[0] << "\n";

	return 1;
}

teuk_blc::teuk_blc(const hbl_parameters params): _params(params) {}
void teuk_blc::operator()(const state_type &y, state_type &dydr, const double r) const
{
	double FR = teuk_potential_FR(r, _params);
	double GR = teuk_potential_GR(r, _params);
	double GI = teuk_potential_GI(r, _params);

	dydr[0] = y[2];
	dydr[1] = y[3];
	dydr[2] = -(FR*y[2] + GR*y[0] - GI*y[1]);
	dydr[3] = -(FR*y[3] + GR*y[1] + GI*y[0]);
}

int teuk_blc_gsl(double r, const double y[], double dydr[], void* parameters){
    hbl_parameters params = *(hbl_parameters *)parameters;

	double FR = teuk_potential_FR(r, params);
	double GR = teuk_potential_GR(r, params);
	double GI = teuk_potential_GI(r, params);

	dydr[0] = y[2];
	dydr[1] = y[3];
	dydr[2] = -(FR*y[2] + GR*y[0] - GI*y[1]);
	dydr[3] = -(FR*y[3] + GR*y[1] + GI*y[0]);

  return GSL_SUCCESS;
}

double teuk_potential_FR(const double r, hbl_parameters params){
	return 2.*(params.s + 1.)*(r - 1.)/(pow(r, 2) - 2.*r + pow(params.a, 2));
}

double teuk_potential_GR(const double r, hbl_parameters params){
	return pow(((r*r + params.a*params.a)*params.om - params.m*params.a)/(pow(r, 2) - 2.*r + pow(params.a, 2)), 2) - params.la/(pow(r, 2) - 2.*r + pow(params.a, 2));
}

double teuk_potential_GI(const double r, hbl_parameters params){
	return -2.*params.s*(r - 1.)*((r*r + params.a*params.a)*params.om - params.m*params.a)/pow(pow(r, 2) - 2.*r + pow(params.a, 2), 2) + 4*params.s*params.om*r/(pow(r, 2) - 2.*r + pow(params.a, 2));
}


double teuk_delta(const double &r, const hbl_parameters &params){
	return pow(r, 2) - 2.*r + pow(params.a, 2);
}

double teuk_deltaP(const double &r, const hbl_parameters &){
	return 2.*(r - 1.);
}

double teuk_deltaPP(const double &, const hbl_parameters &){
	return 2.;
}

double teuk_K(const double &r, const hbl_parameters &params){
	return (pow(r, 2) + pow(params.a, 2))*params.om - params.a*params.m;
}

double teuk_KP(const double &r, const hbl_parameters &params){
	return 2.*r*params.om;
}

double teuk_KPP(const double &, const hbl_parameters &params){
	return 2.*params.om;
}

double teuk_varpi(const double &r, const hbl_parameters &params){
	return pow(r, 2) + pow(params.a, 2);
}

double teuk_varpiP(const double &r, const hbl_parameters &){
	return 2.*r;
}

//*************************************************************//
// 			Generalized Sasaki-Nakamura transformation 		   //
//*************************************************************//

int teuk_in_GSN_integrate(ComplexVector &R_return, ComplexVector &Rp_return, RadialTeukolsky &teuk, const Vector &r){
	hbl_parameters params = {
		.a = teuk.getBlackHoleSpin(),
		.s = teuk.getSpinWeight(),
		.m = teuk.getAzimuthalModeNumber(),
		.om = teuk.getModeFrequency(),
		.la = teuk.getSpinWeightedSpheroidalEigenvalue(),
		.H = 0
	};

	Complex R = teuk.getBoundarySolution(In).getValue();
	Complex RP = teuk.getBoundaryDerivative(In).getValue();
	double r0 = teuk.getBoundaryPoint(In);

	Complex Psi = teuk_R_to_GSN_X(R, RP, r0, params);
	Complex dPsi = teuk_RP_to_GSN_dX(R, RP, r0, params);
	Complex norm = 1.;

	// std::cout << std::setprecision(15);
	// std::cout << "(RADIALSOLVER) Rin(r = "<< r0 << ") = " << R << "\n";
	// std::cout << "(RADIALSOLVER) Rin'(r = "<< r0 << ") = " << RP << "\n";
	// std::cout << "(RADIALSOLVER) Psi(r = "<< r0 << ") = " << Psi << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< r0 << ") = " << dPsi << "\n";

	// std::cout << "(RADIALSOLVER) Reverse transformation returns with accuracy " << std::abs(1. - R/GSN_X_to_teuk_R(Psi, dPsi, r0, params)) << "\n";
	// std::cout << "(RADIALSOLVER) Reverse transformation returns with accuracy " << std::abs(1. - RP/GSN_dX_to_teuk_RP(Psi, dPsi, r0, params)) << "\n";

	double psiR, psiI, dpsiR, dpsiI;
	psiR = std::real(Psi/norm);
	psiI = std::imag(Psi/norm);
	dpsiR = std::real(dPsi/norm);
	dpsiI = std::imag(dPsi/norm);

	int stepNum = r.size();

	state_type Psi0 {psiR, psiI, dpsiR, dpsiI};
	ComplexVector PsiVec(stepNum);
	ComplexVector dPsiVec(stepNum);
	if(BOOST_INTEGRATE){
		teuk_GSN sys(params);
		teuk_integrate_boost(PsiVec, dPsiVec, sys, Psi0, r0, r);
	}else{
		teuk_integrate_gsl(PsiVec, dPsiVec, &teuk_GSN_gsl, Psi0, r0, r, &params);
	}
	// Complex PsiTemp = norm*PsiVec[stepNum-1];
	// Complex dPsiTemp = norm*dPsiVec[stepNum-1];
	//
	// Complex R_f = GSN_X_to_teuk_R(PsiTemp, dPsiTemp, r[stepNum - 1], params);
	// Complex RP_f = GSN_dX_to_teuk_RP(PsiTemp, dPsiTemp, r[stepNum - 1], params);

	// std::cout << "(RADIALSOLVER) Psi(r = "<< r[stepNum - 1] << ") = " << PsiTemp << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< r[stepNum - 1] << ") = " << dPsiTemp << "\n";
	// std::cout << "(RADIALSOLVER) Rin(r = "<< r[stepNum - 1] << ") = " << R_f << "\n";
	// std::cout << "(RADIALSOLVER) Rin'(r = "<< r[stepNum - 1] << ") = " << RP_f << "\n";
	// std::cout << std::setprecision(6);

	for(int i = 0; i < stepNum; i++){
		R_return[i] = GSN_X_to_teuk_R(norm*PsiVec[i], norm*dPsiVec[i], r[i], params);
		Rp_return[i] = GSN_dX_to_teuk_RP(norm*PsiVec[i], norm*dPsiVec[i], r[i], params);
	}

	return 0;
}

int teuk_up_GSN_integrate(ComplexVector &R_return, ComplexVector &Rp_return, RadialTeukolsky &teuk, const Vector &r){
	hbl_parameters params = {
		.a = teuk.getBlackHoleSpin(),
		.s = teuk.getSpinWeight(),
		.m = teuk.getAzimuthalModeNumber(),
		.om = teuk.getModeFrequency(),
		.la = teuk.getSpinWeightedSpheroidalEigenvalue(),
		.H = 0
	};

	// Complex R = teuk.getBoundarySolution(Up).getValue();
	// Complex RP = teuk.getBoundaryDerivative(Up).getValue();
	double r0 = teuk.getBoundaryPoint(Up);
	if(r0 > 50.*r.back()){
		r0 = 50.*r.back();
	}

	// Complex Psi = teuk_R_to_GSN_X(R, RP, r0, params);
	// Complex dPsi = teuk_RP_to_GSN_dX(R, RP, r0, params);
	Complex Psi = gsn_up_asymptotic_infinity(r0, params).getValue();
	Complex dPsi = gsn_up_derivative_asymptotic_infinity(r0, params).getValue();
	Complex norm = 1.;

	while((isnan(std::abs(Psi)) || isnan(std::abs(dPsi))) && r0 < 1.e6*r.back()){
		r0 *= 1.5;
		Psi = gsn_up_asymptotic_infinity(r0, params).getValue();
		dPsi = gsn_up_derivative_asymptotic_infinity(r0, params).getValue();
		norm = 1.;
	}


	// std::cout << std::setprecision(15);
	// std::cout << "(RADIALSOLVER) Rup(r = "<< r0 << ") = " << GSN_X_to_teuk_R(norm*Psi, norm*dPsi, r0, params) << "\n";
	// std::cout << "(RADIALSOLVER) Rup'(r = "<< r0 << ") = " << GSN_dX_to_teuk_RP(norm*Psi, norm*dPsi, r0, params) << "\n";
	// std::cout << "(RADIALSOLVER) Psi(r = "<< r0 << ") = " << Psi << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< r0 << ") = " << dPsi << "\n";

	double psiR, psiI, dpsiR, dpsiI;
	psiR = std::real(Psi/norm);
	psiI = std::imag(Psi/norm);
	dpsiR = std::real(dPsi/norm);
	dpsiI = std::imag(dPsi/norm);

	int stepNum = r.size();

	state_type Psi0 {psiR, psiI, dpsiR, dpsiI};
	ComplexVector PsiVec(stepNum);
	ComplexVector dPsiVec(stepNum);
	Vector rReverse(stepNum);
	for(int i = 0; i < stepNum; i++){
		rReverse[i] = r[stepNum - 1 - i];
	}

	if(BOOST_INTEGRATE){
		teuk_GSN sys(params);
		teuk_integrate_boost(PsiVec, dPsiVec, sys, Psi0, r0, rReverse);
	}else{
		teuk_integrate_gsl(PsiVec, dPsiVec, &teuk_GSN_gsl, Psi0, r0, rReverse, &params);
	}
	// Complex PsiTemp = norm*PsiVec[stepNum-1];
	// Complex dPsiTemp = norm*dPsiVec[stepNum-1];
	//
	// Complex R_f = GSN_X_to_teuk_R(PsiTemp, dPsiTemp, rReverse[stepNum - 1], params);
	// Complex RP_f = GSN_dX_to_teuk_RP(PsiTemp, dPsiTemp, rReverse[stepNum - 1], params);

	// std::cout << "(RADIALSOLVER) Psi(r = "<< rReverse[stepNum - 1] << ") = " << PsiTemp << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< rReverse[stepNum - 1] << ") = " << dPsiTemp << "\n";
	// std::cout << "(RADIALSOLVER) Rup(r = "<< rReverse[stepNum - 1] << ") = " << R_f << "\n";
	// std::cout << "(RADIALSOLVER) Rup'(r = "<< rReverse[stepNum - 1] << ") = " << RP_f << "\n";
	// std::cout << std::setprecision(6);

	for(int i = 0; i < stepNum; i++){
		R_return[i] = GSN_X_to_teuk_R(norm*PsiVec[stepNum - 1 - i], norm*dPsiVec[stepNum - 1 - i], rReverse[stepNum - 1 - i], params);
		Rp_return[i] = GSN_dX_to_teuk_RP(norm*PsiVec[stepNum - 1 - i], norm*dPsiVec[stepNum - 1 - i], rReverse[stepNum - 1 - i], params);
	}

	return 0;
}

teuk_GSN::teuk_GSN(const hbl_parameters params): _params(params) {}
void teuk_GSN::operator()(const state_type &y, state_type &dydr, const double r) const
{
	Complex F = GSN_F_potential_func(r, _params);
	Complex U = GSN_U_potential_func(r, _params);

	dydr[0] = y[2];
	dydr[1] = y[3];
	dydr[2] = std::real(F)*y[2] - std::imag(F)*y[3] + std::real(U)*y[0] - std::imag(U)*y[1];
	dydr[3] = std::real(F)*y[3] + std::imag(F)*y[2] + std::real(U)*y[1] + std::imag(U)*y[0];
}

int teuk_GSN_gsl(double r, const double y[], double dydr[], void* parameters){
    hbl_parameters params = *(hbl_parameters *)parameters;

		Complex F = GSN_F_potential_func(r, params);
		Complex U = GSN_U_potential_func(r, params);

		dydr[0] = y[2];
		dydr[1] = y[3];
		dydr[2] = std::real(F)*y[2] - std::imag(F)*y[3] + std::real(U)*y[0] - std::imag(U)*y[1];
		dydr[3] = std::real(F)*y[3] + std::imag(F)*y[2] + std::real(U)*y[1] + std::imag(U)*y[0];

    return GSL_SUCCESS;
}

// ComplexVector GSN_X_to_teuk_R(const ComplexVector X, const ComplexVector dX, const Vector &r, const hbl_parameters &params){
// 	double a = params.a;
// 	double om = params.om;
// 	double la = params.la;
// 	double m = double(params.m);
// 	ComplexVector Rteuk(r.size());
//
// 	for(int i = 0; i < r.size(); i++){
// 		Rteuk[i] = (2.*I*dX[i]*r[i]*(pow(a, 2) + (-2. + r[i])*r[i])*(a*a + r[i]*r[i])*(-a*m*r[i] - I*(-3. + r[i])*r[i] +
// 	     om*pow(r[i], 3) + a*a*(-2.*I + om*r[i])) + (2.*I*om*(-3. + r[i])*pow(r[i], 6) -
// 	     2.*pow(om, 2)*pow(r[i], 8) + (-2. + r[i])*pow(r[i], 4)*(-6. + (2. + la)*r[i]) +
// 	     4.*pow(a, 5)*m*r[i]*(-I + om*r[i]) + 2.*a*m*pow(r[i], 4)*(I + 2.*om*r[i]*r[i]) +
// 	     2.*pow(a, 3)*m*pow(r[i], 2)*(3.*I - 2.*I*r[i] + 4.*om*r[i]*r[i]) +
// 	     pow(a, 6)*(10. - 2.*om*r[i]*(-2.*I + om*r[i])) + pow(a, 4)*r[i]*(-32. +
// 	        r[i]*(20. + la - 2.*m*m - 2.*om*(3.*I + r[i]*(-5.*I + 3.*om*r[i])))) -
// 	     2.*a*a*r[i]*r[i]*(-12. + r[i]*(19. + la + r[i]*(-6. - la + m*m +
// 				 om*(6.*I + r[i]*(-4.*I + 3.*om*r[i]))))))*X[i])/pow(r[i], 2)/pow(r[i]*r[i] + a*a, 1.5);
// 		Rteuk[i] /= GSN_eta(r[i], params);
// 	}
//
// 	return Rteuk;
// }

Complex teuk_R_to_GSN_X(const Complex R, Complex Rp, const double &r, const hbl_parameters &params){
	return GSN_chi_to_GSN_X(teuk_R_to_GSN_chi(R, Rp, r, params), r, params);
}
Complex teuk_RP_to_GSN_dX(const Complex R, Complex Rp, const double &r, const hbl_parameters &params){
	return GSN_dchi_to_GSN_dX(teuk_R_to_GSN_chi(R, Rp, r, params), teuk_RP_to_GSN_dchi(R, Rp, r, params), r, params);
}
Complex GSN_X_to_teuk_R(const Complex X, Complex dX, const double &r, const hbl_parameters &params){
	return GSN_chi_to_teuk_R(GSN_X_to_GSN_chi(X, dX, r, params), GSN_dX_to_GSN_dchi(X, dX, r, params), r, params);
}
Complex GSN_dX_to_teuk_RP(const Complex X, Complex dX, const double &r, const hbl_parameters &params){
	return GSN_dchi_to_teuk_RP(GSN_X_to_GSN_chi(X, dX, r, params), GSN_dX_to_GSN_dchi(X, dX, r, params), r, params);
}

Complex teuk_R_to_GSN_chi(const Complex R, const Complex Rp, const double &r, const hbl_parameters &params){
	return GSN_alpha(r, params)*R + GSN_beta(r, params)*pow(teuk_delta(r, params), params.s + 1)*Rp;
}
Complex teuk_RP_to_GSN_dchi(const Complex R, const Complex Rp, const double &r, const hbl_parameters &params){
	return ((GSN_alphaP(r, params) + GSN_beta(r, params)*pow(teuk_delta(r, params), params.s)*GSN_V(r, params))*R
		+ (GSN_alpha(r, params) + GSN_betaP(r, params)*pow(teuk_delta(r, params), params.s + 1))*Rp);
}
Complex GSN_chi_to_teuk_R(const Complex chi, const Complex dchi, const double &r, const hbl_parameters &params){
	return ((GSN_alpha(r, params) + GSN_betaP(r, params)*pow(teuk_delta(r, params), params.s + 1))*chi
		- (GSN_beta(r, params)*pow(teuk_delta(r, params), params.s + 1))*dchi)/GSN_eta(r, params);
}
Complex GSN_dchi_to_teuk_RP(const Complex chi, const Complex dchi, const double &r, const hbl_parameters &params){
	return -((GSN_alphaP(r, params) + GSN_beta(r, params)*pow(teuk_delta(r, params), params.s)*GSN_V(r, params))*chi
		- (GSN_alpha(r, params))*dchi)/GSN_eta(r, params);
}

Complex GSN_X_to_GSN_chi(const Complex X, const double &r, const hbl_parameters &params){
	return X/sqrt(teuk_varpi(r, params)*pow(teuk_delta(r, params), params.s));
}
Complex GSN_X_to_GSN_chi(const Complex X, const Complex, const double &r, const hbl_parameters &params){
	return GSN_X_to_GSN_chi(X, r, params);
}
Complex GSN_dX_to_GSN_dchi(const Complex X, const Complex dX, const double &r, const hbl_parameters &params){
	double varpi = teuk_varpi(r, params);
	double varpiP = teuk_varpiP(r, params);
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return (dX/sqrt(varpi*pow(delta, params.s)) - 0.5*(varpiP*pow(delta, params.s) + params.s*varpi*deltaP*pow(delta, params.s - 1))*X/pow(varpi*pow(delta, params.s), 1.5));
}
Complex GSN_chi_to_GSN_X(const Complex chi, const double &r, const hbl_parameters &params){
	return chi*sqrt(teuk_varpi(r, params)*pow(teuk_delta(r, params), params.s));
}
Complex GSN_chi_to_GSN_X(const Complex chi, const Complex, const double &r, const hbl_parameters &params){
	return GSN_chi_to_GSN_X(chi, r, params);
}
Complex GSN_dchi_to_GSN_dX(const Complex chi, const Complex dchi, const double &r, const hbl_parameters &params){
	double varpi = teuk_varpi(r, params);
	double varpiP = teuk_varpiP(r, params);
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return (dchi*sqrt(varpi*pow(delta, params.s)) + 0.5*(varpiP*pow(delta, params.s) + params.s*varpi*deltaP*pow(delta, params.s - 1))*chi/sqrt(varpi*pow(delta, params.s)));
}

Result GSN_chi_to_GSN_X(const Result chi, const double &r, const hbl_parameters &params){
	return chi*sqrt(teuk_varpi(r, params)*pow(teuk_delta(r, params), params.s));
}
Result GSN_chi_to_GSN_X(const Result chi, const Result, const double &r, const hbl_parameters &params){
	return GSN_chi_to_GSN_X(chi, r, params);
}
Result GSN_dchi_to_GSN_dX(const Result chi, const Result dchi, const double &r, const hbl_parameters &params){
	double varpi = teuk_varpi(r, params);
	double varpiP = teuk_varpiP(r, params);
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return (dchi*sqrt(varpi*pow(delta, params.s)) + 0.5*(varpiP*pow(delta, params.s) + params.s*varpi*deltaP*pow(delta, params.s - 1))*chi/sqrt(varpi*pow(delta, params.s)));
}

Complex GSN_eta(const double &r, const hbl_parameters &params){
	Complex c = params.la*(params.la + 2.) - 12.*params.a*params.om*(params.a*params.om - params.m);
	c += (12.*pow(params.a, 2)*(1. - 2.*pow(params.a*params.om - params.m, 2)))*pow(r, -2);
	c += -24.*pow(params.a, 2)*pow(r, -3);
	c += 12.*pow(params.a, 4)*pow(r, -4);

	c += -12.*I*params.om;
	c += 8.*I*params.a*(3*params.a*params.om - params.la*(params.a*params.om - params.m))*pow(r, -1);
	c += -24.*I*params.a*(params.a*params.om - params.m)*pow(r, -2);
	c += 24.*I*pow(params.a, 3)*(params.a*params.om - params.m)*pow(r, -3);

	return c;
}

Complex GSN_etaP(const double &r, const hbl_parameters &params){
	Complex c = 12.*pow(params.a, 2)*(1. - 2.*pow(params.a*params.om - params.m, 2))*2.*pow(r, -2);
	c += -24.*pow(params.a, 2)*3.*pow(r, -3);
	c += 12.*pow(params.a, 4)*4.*pow(r, -4);

	c += 8.*I*params.a*(3*params.a*params.om - params.la*(params.a*params.om - params.m))*pow(r, -1);
	c += -24.*I*params.a*(params.a*params.om - params.m)*2.*pow(r, -2);
	c += 24.*I*pow(params.a, 3)*(params.a*params.om - params.m)*3.*pow(r, -3);

	return -c/r;
}

// Complex GSN_U_potential_func(const double &r, const hbl_parameters &params){
// 	return pow(teuk_varpi(r, params)/teuk_delta(r, params), 2)*GSN_U_star_potential(r, params);
// }
//
// Complex GSN_F_potential_func(const double &r, const hbl_parameters &params){
// 	return -(2.*(pow(r, 2) - pow(params.a, 2))
// 		- pow(teuk_varpi(r, params), 2)*GSN_F_star_potential(r, params))/(teuk_delta(r, params)*teuk_varpi(r, params));
// }

Complex GSN_U_potential_func(const double &r, const hbl_parameters &params){
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	Complex beta = GSN_beta(r, delta, deltaP, params);
	Complex betaP = GSN_betaP(r, delta, deltaP, params);
	Complex betaPP = GSN_betaPP(r, delta, deltaP, params);

	return pow(teuk_varpi(r, params)/delta, 2)*GSN_U_star_potential(r, delta, deltaP, beta, betaP, betaPP, params);
}

Complex GSN_F_potential_func(const double &r, const hbl_parameters &params){
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return -(2.*(pow(r, 2) - pow(params.a, 2))
		- pow(teuk_varpi(r, params), 2)*GSN_F_star_potential(r, delta, deltaP, params))/(delta*teuk_varpi(r, params));
}

Complex GSN_F_star_potential(const double &r, const double &delta, const double &, const hbl_parameters &params){
	return GSN_etaP(r, params)/GSN_eta(r, params)*delta/(r*r + params.a*params.a);
}

Complex GSN_U_star_potential(const double &r, const double &delta, const double &deltaP, const Complex &beta, const Complex &betaP, const Complex &betaPP, const hbl_parameters &params){
	return (GSN_U1(r, delta, deltaP, beta, betaP, betaPP, params)*delta/pow((r*r + params.a*params.a), 2)
		+ pow(GSN_G(r, delta, deltaP, params), 2) + GSN_GP(r, delta, deltaP, params)*delta/(r*r + params.a*params.a)
			 - GSN_G(r, delta, deltaP, params)*GSN_F_star_potential(r, delta, deltaP, params));
}

Complex GSN_G(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params){
	return 0.5*params.s*deltaP/(r*r + params.a*params.a) + r*delta/pow(r*r + params.a*params.a, 2);
}

Complex GSN_GP(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params){
	double varpi = (r*r + params.a*params.a);
	double varpiP = 2.*r;

	return (params.s/varpi - 0.5*params.s*deltaP*varpiP/pow(varpi, 2) + delta/pow(varpi, 2) + r*deltaP/pow(varpi, 2)
			- 2.*r*delta*varpiP/pow(varpi, 3));
}

Complex GSN_U1(const double &r, const double &delta, const double &deltaP, const Complex &beta, const Complex &betaP, const Complex &betaPP, const hbl_parameters &params){
	return GSN_V(r, params) + pow(delta, 2)/beta*(2.*GSN_alphaP(r, delta, deltaP, beta, betaP, params)
		+ betaPP/delta - betaP*deltaP/pow(delta, 2)
			- GSN_etaP(r, params)/GSN_eta(r, params)*(GSN_alpha(r, delta, deltaP, beta, params) + betaP/delta));
}

Complex GSN_alpha(const double &r, const double &delta, const double &, const Complex &beta, const hbl_parameters &params){
	return -I*teuk_K(r, params)*beta/pow(delta, 2) + 3.*I*teuk_KP(r, params) + params.la + 6.*delta/pow(r, 2);
}

Complex GSN_alphaP(const double &r, const double &delta, const double &deltaP, const Complex &beta, const Complex &betaP, const hbl_parameters &params){
	return (-I*teuk_KP(r, params)*beta/pow(delta, 2)
		-I*teuk_K(r, params)*betaP/pow(delta, 2)
			+ 2.*I*teuk_K(r, params)*beta*deltaP/pow(delta, 3)
				+ 6.*I*params.om + 6.*deltaP/pow(r, 2) - 12.*delta/pow(r, 3));
}

Complex GSN_beta(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params){
	return delta*(deltaP - 2.*I*teuk_K(r, params) - 4*delta/r);
}

Complex GSN_betaP(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params){
	return (deltaP*(deltaP - 2.*I*teuk_K(r, params) - 4*delta/r)
		+ delta*(2. - 2.*I*teuk_KP(r, params) - 4*deltaP/r + 4*delta/pow(r, 2)));
}

Complex GSN_betaPP(const double &r, const double &delta, const double &deltaP, const hbl_parameters &params){
	return (2.*(deltaP - 2.*I*teuk_K(r, params) - 4*delta/r)
		+ 2.*deltaP*(2. - 2.*I*teuk_KP(r, params) - 4*deltaP/r + 4*delta/pow(r, 2))
			+ delta*(-4.*I*params.om - 8./r + 8*deltaP/pow(r, 2) - 8*delta/pow(r, 3)));
}

Complex GSN_alpha(const double &r, const hbl_parameters &params){
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return GSN_alpha(r, delta, deltaP, GSN_beta(r, delta, deltaP, params), params);
}

Complex GSN_alphaP(const double &r, const hbl_parameters &params){
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return GSN_alphaP(r, delta, deltaP, GSN_beta(r, delta, deltaP, params), GSN_betaP(r, delta, deltaP, params), params);
}

Complex GSN_beta(const double &r, const hbl_parameters &params){
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return GSN_beta(r, delta, deltaP, params);
}
Complex GSN_betaP(const double &r, const hbl_parameters &params){
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return GSN_betaP(r, delta, deltaP, params);
}
Complex GSN_betaPP(const double &r, const hbl_parameters &params){
	double delta = teuk_delta(r, params);
	double deltaP = teuk_deltaP(r, params);
	return GSN_betaPP(r, delta, deltaP, params);
}

Complex GSN_V(const double &r, const hbl_parameters &params){
	double K = teuk_K(r, params);
	return params.la - 4.*I*Complex(params.s)*params.om*r - (pow(K, 2) - 2.*I*Complex(params.s)*(r - 1.)*K)/(r*r - 2.*r + params.a*params.a);
}


Complex GSN_U_potential(const double &r, const hbl_parameters &params){
	return GSN_U_potential(params.a, params.s, params.m, params.om, params.la, r);
}

Complex GSN_F_potential(const double &r, const hbl_parameters &params){
	return GSN_F_potential(params.a, params.s, params.m, params.om, params.la, r);
}

Complex GSN_U_potential(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &r){
	Complex m = Complex(mTemp);
	return pow((-2. + r)*r + pow(a,2),-2)*pow(pow(a,2) + pow(r,2),-2)*(((-2. + r)*r + pow(a,2))*(-8.*r*pow(a,2) - 1.*pow(a,4) + pow(r,4)) + pow((-2. + r)*pow(a,2) + pow(r,3),2) - 8.*a*((-2. + r)*r + pow(a,2))*pow(r,-1)*(pow(a,2) + pow(r,2))*(3.*m*r*(Complex(0.,3.) - 4.*omega*r)*pow(a,2) - Complex(0.,1.)*m*(6. + lambda*r)*pow(r,2) + a*r*(9. + r*(-3. + Complex(0.,6.)*omega + 6.*pow(m,2)) + Complex(0.,1.)*(-3. + lambda)*omega*pow(r,2)) + pow(a,3)*(-6. - Complex(0.,9.)*omega*r + 6.*pow(omega,2)*pow(r,2)))*((-2. + r)*pow(a,2) + pow(r,3))*pow(24.*m*r*(Complex(0.,1.) - 2.*omega*r)*pow(a,3) - 4.*a*m*pow(r,2)*(Complex(0.,6.) + Complex(0.,2.)*lambda*r + 3.*omega*pow(r,2)) + 12.*pow(a,4)*(-1. - Complex(0.,2.)*omega*r + 2.*pow(omega,2)*pow(r,2)) + 4.*r*pow(a,2)*(6. + r*(-3. + Complex(0.,6.)*omega + 6.*pow(m,2)) + Complex(0.,2.)*(-3. + lambda)*omega*pow(r,2) + 3.*pow(omega,2)*pow(r,3)) - 1.*(2.*lambda - Complex(0.,12.)*omega + pow(lambda,2))*pow(r,4),-1) + ((-2. + r)*r + pow(a,2))*pow(pow(a,2) + pow(r,2),2)*(lambda + Complex(0.,8.)*omega*r - 1.*(-1.*a*m + omega*(pow(a,2) + pow(r,2)))*(-1.*a*m + Complex(0.,4.)*(-1. + r) + omega*(pow(a,2) + pow(r,2)))*pow((-2. + r)*r + pow(a,2),-1) + 4.*pow(r,-2)*(12.*m*omega*(-1. + Complex(0.,3.)*omega*r)*pow(a,7)*pow(r,2) + Complex(0.,2.)*m*pow(a,5)*pow(r,2)*(-15. + r*(15. - 7.*lambda - Complex(0.,48.)*omega + 6.*pow(m,2)) + Complex(0.,2.)*(15. + lambda)*omega*pow(r,2) + 42.*pow(omega,2)*pow(r,3)) + 6.*pow(a,8)*(1. + pow(omega,2)*pow(r,2) - Complex(0.,2.)*pow(omega,3)*pow(r,3)) - Complex(0.,2.)*m*pow(a,3)*pow(r,3)*(-30. + (33. - 22.*lambda - Complex(0.,24.)*omega)*r + (-3. + 12.*lambda - Complex(0.,6.)*omega - Complex(0.,4.)*lambda*omega + 6.*pow(m,2))*pow(r,2) - Complex(0.,3.)*(2. + lambda + Complex(0.,28.)*omega)*omega*pow(r,3) - 36.*pow(omega,2)*pow(r,4)) + 2.*pow(a,4)*pow(r,3)*(-6.*(-4. + lambda + Complex(0.,5.)*omega) + r*(-9. + lambda + Complex(0.,48.)*omega - Complex(0.,22.)*lambda*omega - 1.*pow(lambda,2) + 12.*pow(omega,2)) + 2.*omega*(Complex(0.,-6.) + 12.*omega + lambda*(Complex(0.,6.) + omega))*pow(r,2) + pow(m,2)*(-24. + (27. + lambda)*r - Complex(0.,12.)*omega*pow(r,2)) + 3.*(-7. + 2.*lambda + Complex(0.,14.)*omega)*pow(omega,2)*pow(r,3) - Complex(0.,24.)*pow(omega,3)*pow(r,4)) + r*pow(a,6)*(-12. + 3.*r*(-4. + lambda + Complex(0.,10.)*omega + 2.*pow(m,2)) - 2.*omega*(Complex(0.,18.) - Complex(0.,7.)*lambda + 24.*omega + Complex(0.,18.)*pow(m,2))*pow(r,2) + 2.*(3. + lambda)*pow(omega,2)*pow(r,3) - Complex(0.,48.)*pow(omega,3)*pow(r,4)) + Complex(0.,1.)*a*m*pow(r,4)*(-24. - 24.*(-1. + lambda)*r - 1.*(-18.*lambda + Complex(0.,48.)*omega + pow(lambda,2))*pow(r,2) + lambda*(lambda - Complex(0.,28.)*omega)*pow(r,3) + 2.*(Complex(0.,5.)*lambda - 18.*omega)*omega*pow(r,4) + 12.*pow(omega,2)*pow(r,5)) + pow(a,2)*pow(r,4)*(12.*(-2. + lambda + Complex(0.,2.)*omega) + 6.*r*(2. - Complex(0.,10.)*omega + Complex(0.,4.)*lambda*omega + pow(lambda,2)) + Complex(0.,1.)*(lambda*(Complex(0.,3.) - 18.*omega) + 6.*(1. + Complex(0.,8.)*omega)*omega + (Complex(0.,3.) + omega)*pow(lambda,2))*pow(r,2) + omega*(Complex(0.,12.) + Complex(0.,4.)*lambda + 108.*omega - 28.*lambda*omega + Complex(0.,1.)*pow(lambda,2))*pow(r,3) + pow(m,2)*(24. + 4.*(-9. + lambda)*r - 6.*(lambda - Complex(0.,14.)*omega)*pow(r,2) - Complex(0.,24.)*omega*pow(r,3)) + 2.*(-15. + 5.*lambda + Complex(0.,18.)*omega)*pow(omega,2)*pow(r,4) - Complex(0.,12.)*pow(omega,3)*pow(r,5)) + (Complex(0.,2.)*lambda + 12.*omega + Complex(0.,1.)*pow(lambda,2))*(Complex(0.,3.) - Complex(0.,2.)*r - 3.*omega*pow(r,2) + omega*pow(r,3))*pow(r,6))*pow((-2. + r)*r + pow(a,2),-1)*pow(24.*m*r*(Complex(0.,1.) - 2.*omega*r)*pow(a,3) - 4.*a*m*pow(r,2)*(Complex(0.,6.) + Complex(0.,2.)*lambda*r + 3.*omega*pow(r,2)) + 12.*pow(a,4)*(-1. - Complex(0.,2.)*omega*r + 2.*pow(omega,2)*pow(r,2)) + 4.*r*pow(a,2)*(6. + r*(-3. + Complex(0.,6.)*omega + 6.*pow(m,2)) + Complex(0.,2.)*(-3. + lambda)*omega*pow(r,2) + 3.*pow(omega,2)*pow(r,3)) - 1.*(2.*lambda - Complex(0.,12.)*omega + pow(lambda,2))*pow(r,4),-1)));
}

Complex GSN_F_potential(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &r){
	Complex m = Complex(mTemp);
	return 2.*pow(r,-1)*(12.*m*r*(Complex(0.,-3.) + 4.*omega*r)*pow(a,7) + 4.*m*pow(a,5)*pow(r,2)*(Complex(0.,30.) + Complex(0.,1.)*(-18. + lambda + Complex(0.,36.)*omega)*r + 24.*omega*pow(r,2)) + pow(a,8)*(24. + Complex(0.,36.)*omega*r - 24.*pow(omega,2)*pow(r,2)) +
     4.*m*pow(a,3)*pow(r,3)*(Complex(0.,-18.) - Complex(0.,4.)*(-6. + lambda)*r + Complex(0.,1.)*(-9. + 2.*lambda + Complex(0.,15.)*omega)*pow(r,2) + 12.*omega*pow(r,3)) -
     4.*r*pow(a,6)*(24. + 3.*r*(-5. + Complex(0.,10.)*omega + 2.*pow(m,2)) + Complex(0.,1.)*(-21. + lambda + Complex(0.,18.)*omega)*omega*pow(r,2) + 12.*pow(omega,2)*pow(r,3)) +
     pow(a,2)*(48. + r*(-48. - 2.*lambda + Complex(0.,36.)*omega - 1.*pow(lambda,2) + 24.*pow(m,2)) - 12.*(-1. + Complex(0.,2.)*omega + 2.*pow(m,2))*pow(r,2) - 4.*omega*(Complex(0.,-3.) + Complex(0.,1.)*lambda + 3.*omega)*pow(r,3))*pow(r,4) -
     4.*pow(a,4)*pow(r,2)*(-24. - 18.*r*(-2. + Complex(0.,1.)*omega + pow(m,2)) + 4.*(-3. - Complex(0.,1.)*(-9. + lambda)*omega + 3.*pow(m,2))*pow(r,2) + Complex(0.,1.)*(-15. + 2.*lambda + Complex(0.,9.)*omega)*omega*pow(r,3) + 6.*pow(omega,2)*pow(r,4)) +
     4.*a*m*(Complex(0.,-6.) + Complex(0.,6.)*r + (Complex(0.,1.)*lambda + 3.*omega)*pow(r,2))*pow(r,5) + (2.*lambda - Complex(0.,12.)*omega + pow(lambda,2))*pow(r,7))*pow((-2. + r)*r + pow(a,2),-1)*pow(pow(a,2) + pow(r,2),-1)*
   pow(24.*m*r*(Complex(0.,1.) - 2.*omega*r)*pow(a,3) - 4.*a*m*pow(r,2)*(Complex(0.,6.) + Complex(0.,2.)*lambda*r + 3.*omega*pow(r,2)) + 12.*pow(a,4)*(-1. - Complex(0.,2.)*omega*r + 2.*pow(omega,2)*pow(r,2)) +
     4.*r*pow(a,2)*(6. + r*(-3. + Complex(0.,6.)*omega + 6.*pow(m,2)) + Complex(0.,2.)*(-3. + lambda)*omega*pow(r,2) + 3.*pow(omega,2)*pow(r,3)) - 1.*(2.*lambda - Complex(0.,12.)*omega + pow(lambda,2))*pow(r,4),-1);
}

//*************************************************************//
// 				Hyperboloidal slicing transformation 		   //
//*************************************************************//

int teuk_in_HBL_integrate(ComplexVector &R_return, ComplexVector &Rp_return, RadialTeukolsky &teuk, const Vector &r){
	hbl_parameters params = {
		.a = teuk.getBlackHoleSpin(),
		.s = teuk.getSpinWeight(),
		.m = teuk.getAzimuthalModeNumber(),
		.om = teuk.getModeFrequency(),
		.la = teuk.getSpinWeightedSpheroidalEigenvalue(),
		.H = -1
	};

	Complex R = teuk.getBoundarySolution(In).getValue();
	Complex RP = teuk.getBoundaryDerivative(In).getValue();
	double r0 = teuk.getBoundaryPoint(In);
	// if(teuk.getSpinWeight() > 0){
	// 	std::cout << "(RADIALSOLVER) Transform to s = " << - params.s << " equation \n";
	// 	teuk_spin_minus_transform_derivative(RP, R, r0, params);
	// 	teuk_spin_minus_transform_solution(R, r0, params);
	// 	teuk_spin_minus_transform_parameters(params);
	// }

	Complex Psi = teuk_R_to_hbl_Psi(R, r0, params);
	Complex dPsi = teuk_RP_to_hbl_dPsi(RP, R, r0, params);
	Complex norm = 1;

	double psiR, psiI, dpsiR, dpsiI;
	psiR = std::real(Psi/norm);
	psiI = std::imag(Psi/norm);
	dpsiR = std::real(dPsi/norm);
	dpsiI = std::imag(dPsi/norm);

	int stepNum = r.size();

	// std::cout << "(RADIALSOLVER) Rin(r = "<< r0 << ") = " << R << "\n";
	// std::cout << "(RADIALSOLVER) Rin'(r = "<< r0 << ") = " << RP << "\n";
	// std::cout << "(RADIALSOLVER) Psi(r = "<< r0 << ") = " << norm*(psiR + I*psiI) << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< r0 << ") = " << norm*(dpsiR + I*dpsiI) << "\n";

	state_type Psi0 {psiR, psiI, dpsiR, dpsiI};
	ComplexVector PsiVec(stepNum);
	ComplexVector dPsiVec(stepNum);
	Vector rReverse = r;
	if(r0 >= r.back()){
		std::reverse(rReverse.begin(),rReverse.end());
	}

	if(BOOST_INTEGRATE){
		teuk_hbl sys(params);
		teuk_integrate_boost(PsiVec, dPsiVec, sys, Psi0, r0, rReverse);
	}else{
		teuk_integrate_gsl(PsiVec, dPsiVec, &teuk_hbl_gsl, Psi0, r0, rReverse, &params);
	}
	// Complex PsiTemp = norm*PsiVec[stepNum-1];
	// Complex dPsiTemp = norm*dPsiVec[stepNum-1];
	//
	// Complex R_f = hbl_Psi_to_teuk_R(PsiTemp, rReverse[stepNum - 1], params);
	// Complex RP_f = hbl_dPsi_to_teuk_RP(dPsiTemp, PsiTemp, rReverse[stepNum - 1], params);
	//
	// std::cout << "(RADIALSOLVER) Psi(r = "<< rReverse[stepNum - 1] << ") = " << PsiTemp << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< rReverse[stepNum - 1] << ") = " << dPsiTemp << "\n";
	// std::cout << "(RADIALSOLVER) Rin(r = "<< rReverse[stepNum - 1] << ") = " << R_f << "\n";
	// std::cout << "(RADIALSOLVER) Rin'(r = "<< rReverse[stepNum - 1] << ") = " << RP_f << "\n";

	for(int i = 0; i < stepNum; i++){
		R_return[i] = hbl_Psi_to_teuk_R(norm*PsiVec[i], r[i], params);
		Rp_return[i] = hbl_dPsi_to_teuk_RP(norm*dPsiVec[i], norm*PsiVec[i], r[i], params);
	}

	if(r0 >= r.back()){
		std::reverse(R_return.begin(),R_return.end());
		std::reverse(Rp_return.begin(),Rp_return.end());
	}

	// if(teuk.getSpinWeight() > 0){
	// 	for(int i = 0; i < r.size(); i++){
	// 		teuk_spin_minus_transform_derivative(Rp_return[i], R_return[i], r[i], params);
	// 		teuk_spin_minus_transform_solution(R_return[i], r[i], params);
	// 	}
	// }

	return 0;
}

int teuk_up_HBL_integrate(ComplexVector &R_return, ComplexVector &Rp_return, RadialTeukolsky &teuk, const Vector &r){
	hbl_parameters params = {
		.a = teuk.getBlackHoleSpin(),
		.s = teuk.getSpinWeight(),
		.m = teuk.getAzimuthalModeNumber(),
		.om = teuk.getModeFrequency(),
		.la = teuk.getSpinWeightedSpheroidalEigenvalue(),
		.H = 1
	};

	Complex R = teuk.getBoundarySolution(Up).getValue();
	Complex RP = teuk.getBoundaryDerivative(Up).getValue();
	double r0 = teuk.getBoundaryPoint(Up);

	Complex Psi = teuk_R_to_hbl_Psi(R, r0, params);
	Complex dPsi = teuk_RP_to_hbl_dPsi(RP, R, r0, params);
	Complex norm = 1.;

	// std::cout << "(RADIALSOLVER) Rup(r = "<< r0 << ") = " << R << "\n";
	// std::cout << "(RADIALSOLVER) Rup'(r = "<< r0 << ") = " << RP << "\n";
	// std::cout << "(RADIALSOLVER) Psi(r = "<< r0 << ") = " << Psi << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< r0 << ") = " << dPsi << "\n";

	double psiR, psiI, dpsiR, dpsiI;
	psiR = std::real(Psi/norm);
	psiI = std::imag(Psi/norm);
	dpsiR = std::real(dPsi/norm);
	dpsiI = std::imag(dPsi/norm);

	int stepNum = r.size();

	state_type Psi0 {psiR, psiI, dpsiR, dpsiI};
	ComplexVector PsiVec(stepNum);
	ComplexVector dPsiVec(stepNum);
	Vector rReverse = r;
	if(r0 >= r.back()){
		std::reverse(rReverse.begin(),rReverse.end());
	}

	if(BOOST_INTEGRATE){
		teuk_hbl sys(params);
		teuk_integrate_boost(PsiVec, dPsiVec, sys, Psi0, r0, rReverse);
	}else{
		teuk_integrate_gsl(PsiVec, dPsiVec, &teuk_hbl_gsl, Psi0, r0, rReverse, &params);
	// 	int success = teuk_integrate_gsl(PsiVec, dPsiVec, &teuk_hbl_gsl, &jacobian, Psi0, r0, rReverse, &params);
	}

	// std::cout << "(RADIALSOLVER) Psi(r = "<< rReverse[0] << ") = " << PsiVec[0] << "\n";
	// std::cout << "(RADIALSOLVER) Psi'(r = "<< rReverse[0] << ") = " << dPsiVec[0] << "\n";

	for(int i = 0; i < stepNum; i++){
		R_return[i] = hbl_Psi_to_teuk_R(norm*PsiVec[i], rReverse[i], params);
		Rp_return[i] = hbl_dPsi_to_teuk_RP(norm*dPsiVec[i], norm*PsiVec[i], rReverse[i], params);
	}

	if(r0 >= r.back()){
		std::reverse(R_return.begin(), R_return.end());
		std::reverse(Rp_return.begin(), Rp_return.end());
	}

	return 0;
}


teuk_hbl::teuk_hbl(const hbl_parameters params): _params(params) {}
void teuk_hbl::operator()(const state_type &y, state_type &dydr, const double r) const
{
	double a = _params.a;
	double delta = r*r - 2.*r + a*a;
	double GR = potential_GR(r, _params);
	double GI = potential_GI(r, _params);
	double UR = potential_UR(r, _params);
	double UI = potential_UI(r, _params);

	//std::cout << "GR("<<r<<") = " << GR << "\n";

	dydr[0] = y[2];
	dydr[1] = y[3];
	dydr[2] = (2.*r*(GR*y[2] - GI*y[3]) - (UR*y[0] - UI*y[1]))/(delta*delta*r*r);
	dydr[3] = (2.*r*(GR*y[3] + GI*y[2]) - (UR*y[1] + UI*y[0]))/(delta*delta*r*r);
}

teuk_hbl_implicit::teuk_hbl_implicit(const hbl_parameters params): _params(params) {}
void teuk_hbl_implicit::operator()(const vector_type &y, vector_type &dydr, const double r) const
{
	double a = _params.a;
	double delta = r*r - 2.*r + a*a;
	double GR = potential_GR(r, _params);
	double GI = potential_GI(r, _params);
	double UR = potential_UR(r, _params);
	double UI = potential_UI(r, _params);

	//std::cout << "GR("<<r<<") = " << GR << "\n";

	dydr[0] = y[2];
	dydr[1] = y[3];
	dydr[2] = (2.*r*(GR*y[2] - GI*y[3]) - (UR*y[0] - UI*y[1]))/(delta*delta*r*r);
	dydr[3] = (2.*r*(GR*y[3] + GI*y[2]) - (UR*y[1] + UI*y[0]))/(delta*delta*r*r);
}

jacobi_hbl_implicit::jacobi_hbl_implicit(const hbl_parameters params): _params(params) {}
void jacobi_hbl_implicit::operator()( const vector_type &y, matrix_type &dfdy, const double r, vector_type &dfdr) const
{
	double a = _params.a;
	double delta = r*r - 2.*r + a*a;
	double GR = potential_GR(r, _params);
	double GI = potential_GI(r, _params);
	double UR = potential_UR(r, _params);
	double UI = potential_UI(r, _params);

	double ddelta = 2.*(r - 1.);
	double dGR = potential_dGR(r, _params);
	double dGI = potential_dGI(r, _params);
	double dUR = potential_dUR(r, _params);
	double dUI = potential_dUI(r, _params);

	dfdy(0, 0) = 0.0; // df[0]/dy[0]
	dfdy(0, 1) = 0.0; // df[0]/dy[1]
	dfdy(0, 2) = 1.0; // df[0]/dy[2]
	dfdy(0, 3) = 0.0; // df[0]/dy[3]

	dfdy(1, 0) = 0.0; // df[1]/dy[0]
	dfdy(1, 1) = 0.0; // df[1]/dy[1]
	dfdy(1, 2) = 0.0; // df[1]/dy[2]
	dfdy(1, 3) = 1.0; // df[1]/dy[3]

	dfdy(2, 0) = -UR/(delta*delta*r*r); // df[2]/dy[0]
	dfdy(2, 1) = UI/(delta*delta*r*r);  // df[2]/dy[1]
	dfdy(2, 2) = 2.*GR/(delta*delta*r);  // df[2]/dy[2]
	dfdy(2, 3) = -2.*GI/(delta*delta*r); // df[2]/dy[3]

	dfdy(3, 0) = -UI/(delta*delta*r*r); // df[3]/dy[0]
	dfdy(3, 1) = -UR/(delta*delta*r*r); // df[3]/dy[1]
	dfdy(3, 2) = 2.*GI/(delta*delta*r);  // df[3]/dy[2]
	dfdy(3, 3) = 2.*GR/(delta*delta*r);  // df[3]/dy[3]

	dfdr[0] = 0.0;
	dfdr[1] = 0.0;
	dfdr[2] = 2.*(UR*y[0] - UI*y[1])/(pow(delta, 2)*pow(r, 3))
		+ 2.*ddelta*(UR*y[0] - UI*y[1])/(pow(delta, 3)*pow(r, 2))
		- (dUR*y[0] - dUI*y[1])/(pow(delta, 2)*pow(r, 2))
		- 2.*(GR*y[2] - GI*y[3])/(pow(delta, 2)*pow(r, 2))
		- 4.*ddelta*(GR*y[2] - GI*y[3])/(pow(delta, 3)*pow(r, 1))
		+ 2.*(dGR*y[2] - dGI*y[3])/(pow(delta, 2)*pow(r, 1));
	dfdr[3] = 2.*(UI*y[0] + UR*y[1])/(pow(delta, 2)*pow(r, 3))
		+ 2.*ddelta*(UI*y[0] + UR*y[1])/(pow(delta, 3)*pow(r, 2))
		- (dUI*y[0] + dUR*y[1])/(pow(delta, 2)*pow(r, 2))
		- 2.*(GI*y[2] + GR*y[3])/(pow(delta, 2)*pow(r, 2))
		- 4.*ddelta*(GI*y[2] + GR*y[3])/(pow(delta, 3)*pow(r, 1))
		+ 2*(dGI*y[2] + dGR*y[3])/(pow(delta, 2)*pow(r, 1));
};

double PhiHBL(double r, double a){
	double kappa = sqrt(1. - a*a);
	return a/(2.*kappa)*log((r - (1. + kappa))/(r - (1. - kappa)));
}

double hHBL(double r, double a){
	return r_to_rstar(r, a);
}

double r_to_rstar(double r, double a){
	double kappa = sqrt(1. - a*a);
	return r + (1. + kappa)*log((r - (1. + kappa)))/kappa - (1. - kappa)*log((r - (1. - kappa)))/kappa;
}

Complex teuk_R_to_hbl_Psi(Complex R, double r, hbl_parameters params){
	double a = params.a, om = params.om;
	Complex s = Complex(params.s), m = Complex(params.m);
	int H = params.H;
	double delta = r*r - 2.*r + a*a;
	double h = H*hHBL(r, a);
	double Phi = PhiHBL(r, a);

	Complex prefactor = r*pow(delta, s)*exp(-I*m*Phi)*exp(-I*om*h);
	Complex Psi = prefactor*R;

	return Psi;
}

Complex teuk_RP_to_hbl_dPsi(Complex RinP, Complex Rin, double r, hbl_parameters params){
	double a = params.a, om = params.om;
	Complex s = Complex(params.s), m = Complex(params.m);
	int H = params.H;
	double delta = r*r - 2.*r + a*a, ddelta = 2.*(r - 1.);
	double h = H*hHBL(r, a), dh = H*(r*r + a*a)/delta;
	double Phi = PhiHBL(r, a), dPhi = a/delta;

	Complex prefactor = r*pow(delta, s)*exp(-I*m*Phi)*exp(-I*om*h);
	Complex dprefactor = (-I*dPhi*m - I*dh*om + 1./r + (s*ddelta)/delta);
	Complex dPsi = prefactor*(dprefactor*Rin + RinP);

	return dPsi;
}

Complex hbl_Psi_to_teuk_R(Complex Psi, double r, hbl_parameters params){
	double a = params.a, om = params.om;
	Complex s = Complex(params.s), m = Complex(params.m);
	int H = params.H;
	double delta = r*r - 2.*r + a*a;
	double h = H*hHBL(r, a);
	double Phi = PhiHBL(r, a);

	Complex prefactor = r*pow(delta, s)*exp(-I*m*Phi)*exp(-I*om*h);
	Complex Rin = Psi/prefactor;

	return Rin;
}

Complex hbl_dPsi_to_teuk_RP(Complex dPsi, Complex Psi, double r, hbl_parameters params){
	double a = params.a, om = params.om;
	Complex s = Complex(params.s), m = Complex(params.m);
	int H = params.H;
	double delta = r*r - 2.*r + a*a, ddelta = 2.*(r - 1.);
	double h = H*hHBL(r, a), dh = H*(r*r + a*a)/delta;
	double Phi = PhiHBL(r, a), dPhi = a/delta;

	Complex prefactor = r*pow(delta, s)*exp(-I*m*Phi)*exp(-I*om*h);
	Complex dprefactor = (-I*dPhi*m - I*dh*om + 1./r + (s*ddelta)/delta);
	Complex RinP = (dPsi - dprefactor*Psi)/prefactor;

	return RinP;
}


int teuk_hbl_gsl(double r, const double y[], double dydr[], void* parameters){
  hbl_parameters params = *(hbl_parameters *)parameters;

	double a = params.a;
	double delta = r*r - 2.*r + a*a;
	double GR = potential_GR(r, params);
	double GI = potential_GI(r, params);
	double UR = potential_UR(r, params);
	double UI = potential_UI(r, params);

	dydr[0] = y[2];
	dydr[1] = y[3];
	dydr[2] = (2.*r*(GR*y[2] - GI*y[3]) - (UR*y[0] - UI*y[1]))/(delta*delta*r*r);
	dydr[3] = (2.*r*(GR*y[3] + GI*y[2]) - (UR*y[1] + UI*y[0]))/(delta*delta*r*r);

  return GSL_SUCCESS;
}

int psi_rstar(double rstar, const double y[], double dydr[], void* parameters){
  hbl_parameters params = *(hbl_parameters *)parameters;

	double a = params.a;
	double r = rstar_to_r(rstar, a);
	double varpi = r*r + a*a;

	double FR = potential_FR(r, params);
	double FI = potential_FI(r, params);
	double UR = potential_UR(r, params);
	double UI = potential_UI(r, params);

	dydr[0] = y[2];
	dydr[1] = y[3];
	dydr[2] = (2*r*(FR*y[2]-FI*y[3]) - (UR*y[0]-UI*y[1]))/(varpi*varpi*r*r);
	dydr[3] = (2*r*(FR*y[3]+FI*y[2]) - (UR*y[1]+UI*y[0]))/(varpi*varpi*r*r);

  return GSL_SUCCESS;
}

int jacobian(double r, const double y[], double *dfdy, double dfdr[], void* parameters){
	hbl_parameters params = *(hbl_parameters *)parameters;

	double a = params.a;
	double delta = r*r - 2.*r + a*a;
	double GR = potential_GR(r, params);
	double GI = potential_GI(r, params);
	double UR = potential_UR(r, params);
	double UI = potential_UI(r, params);

	double ddelta = 2.*(r - 1.);
	double dGR = potential_dGR(r, params);
	double dGI = potential_dGI(r, params);
	double dUR = potential_dUR(r, params);
	double dUI = potential_dUI(r, params);

	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
	gsl_matrix* m = &dfdy_mat.matrix;
	gsl_matrix_set(m, 0, 0, 0.0); // df[0]/dy[0]
	gsl_matrix_set(m, 0, 1, 0.0); // df[0]/dy[1]
	gsl_matrix_set(m, 0, 2, 1.0); // df[0]/dy[2]
	gsl_matrix_set(m, 0, 3, 0.0); // df[0]/dy[3]

	gsl_matrix_set(m, 1, 0, 0.0); // df[1]/dy[0]
	gsl_matrix_set(m, 1, 1, 0.0); // df[1]/dy[1]
	gsl_matrix_set(m, 1, 2, 0.0); // df[1]/dy[2]
	gsl_matrix_set(m, 1, 3, 1.0); // df[1]/dy[3]

	gsl_matrix_set(m, 2, 0, -UR/(delta*delta*r*r)); // df[2]/dy[0]
	gsl_matrix_set(m, 2, 1, UI/(delta*delta*r*r));  // df[2]/dy[1]
	gsl_matrix_set(m, 2, 2, 2.*GR/(delta*delta*r));  // df[2]/dy[2]
	gsl_matrix_set(m, 2, 3, -2.*GI/(delta*delta*r)); // df[2]/dy[3]

	gsl_matrix_set(m, 3, 0, -UI/(delta*delta*r*r)); // df[3]/dy[0]
	gsl_matrix_set(m, 3, 1, -UR/(delta*delta*r*r)); // df[3]/dy[1]
	gsl_matrix_set(m, 3, 2, 2.*GI/(delta*delta*r));  // df[3]/dy[2]
	gsl_matrix_set(m, 3, 3, 2.*GR/(delta*delta*r));  // df[3]/dy[3]

	// df[i]/dr
	dfdr[0] = 0.0;
	dfdr[1] = 0.0;
	dfdr[2] = 2.*(UR*y[0] - UI*y[1])/(pow(delta, 2)*pow(r, 3))
		+ 2.*ddelta*(UR*y[0] - UI*y[1])/(pow(delta, 3)*pow(r, 2))
		- (dUR*y[0] - dUI*y[1])/(pow(delta, 2)*pow(r, 2))
		- 2.*(GR*y[2] - GI*y[3])/(pow(delta, 2)*pow(r, 2))
		- 4.*ddelta*(GR*y[2] - GI*y[3])/(pow(delta, 3)*pow(r, 1))
		+ 2.*(dGR*y[2] - dGI*y[3])/(pow(delta, 2)*pow(r, 1));
	dfdr[3] = 2.*(UI*y[0] + UR*y[1])/(pow(delta, 2)*pow(r, 3))
		+ 2.*ddelta*(UI*y[0] + UR*y[1])/(pow(delta, 3)*pow(r, 2))
		- (dUI*y[0] + dUR*y[1])/(pow(delta, 2)*pow(r, 2))
		- 2.*(GI*y[2] + GR*y[3])/(pow(delta, 2)*pow(r, 2))
		- 4.*ddelta*(GI*y[2] + GR*y[3])/(pow(delta, 3)*pow(r, 1))
		+ 2*(dGI*y[2] + dGR*y[3])/(pow(delta, 2)*pow(r, 1));

	return GSL_SUCCESS;
}

double potential_GR(double r, hbl_parameters params){
	double a = params.a;
	int s = params.s;
	double delta = r*r - 2.*r + a*a;
	double GR;

	GR = delta*(pow(a, 2.) + r*(r*s - (1. + s)));

	return GR;
}

double potential_GI(double r, hbl_parameters params){
	double a = params.a, om = params.om;
	int H = params.H, m = params.m;
	double delta = r*r - 2.*r + a*a;
	double varpi = r*r + a*a;
	double GI;

	GI = -(delta*r*(a*m + H*om*varpi));

	return GI;
}

double potential_FR(double r, hbl_parameters params){
	double a = params.a;
	int s = params.s;
	double delta = r*r - 2.*r + a*a;
	double varpi = r*r + a*a;
	double FR;

	FR = delta*pow(a, 2) + varpi*s*r*(r - 1.);

	return FR;
}

double potential_FI(double r, hbl_parameters params){
	double a = params.a, om = params.om;
	int H = params.H, m = params.m;
	double varpi = r*r + a*a;
	double FI;

	FI = -(varpi*r*(a*m + H*om*varpi));

	return FI;
}

double potential_UR(double r, hbl_parameters params){
	double a = params.a, om = params.om, la = params.la;
	int s = params.s, H = params.H, m = params.m;
	double delta = r*r - 2.*r + a*a;
	double varpi = r*r + a*a;
	double UR;

	UR = delta*(2*pow(a, 2) - la*pow(r, 2) - 2.*r*(1. + s))
		- 2.*a*(1. + H)*m*om*pow(r, 2)*varpi
		+ (1. - pow(H, 2))*pow(om, 2)*pow(r, 2)*pow(varpi, 2);

	return UR;
}

double potential_UI(double r, hbl_parameters params){
	double a = params.a, om = params.om;
	int s = params.s, H = params.H, m = params.m;
	double delta = r*r - 2.*r + a*a;
	double UI;

	UI = -2.*a*delta*(m + a*H*om)*r + 2.*om*pow(r, 2)*(delta*(1. - H)*r
		- (1. + H)*(-pow(a, 2) + pow(r, 2)))*s;

	return UI;
}

double potential_dGR(double r, hbl_parameters params){
	double a = params.a;
	int s = params.s;
	double delta = r*r - 2.*r + a*a;
	double ddelta = 2.*(r - 1.);
	double dGR;

	dGR = -delta*(1. + s - 2.*r*s) + ddelta*(pow(a, 2) - r*(1. + s - r*s));

	return dGR;
}

double potential_dGI(double r, hbl_parameters params){
	double a = params.a, om = params.om;
	int H = params.H, m = params.m;
	double delta = r*r - 2.*r + a*a;
	double ddelta = 2.*(r - 1.);
	double varpi = r*r + a*a;
	double dvarpi = 2.*r;
	double dGI;

	dGI = -delta*(a*m + H*om*varpi) - ddelta*r*(a*m + H*om*varpi) - H*om*r*delta*dvarpi;

	return dGI;
}

double potential_dUR(double r, hbl_parameters params){
	double a = params.a, om = params.om, la = params.la;
	int s = params.s, H = params.H, m = params.m;
	double delta = r*r - 2.*r + a*a;
	double ddelta = 2.*(r - 1.);
	double dUR;

	dUR = -2.*(1. + H)*om*r*(pow(a, 3)*(2.*m - a*(1. - H)*om)
		+ 4.*a*(m - a*(1. - H)*om)*pow(r, 2)
		- 3.*(1. - H)*om*pow(r, 4)) - 2.*delta*(1. + la*r + s)
		+ ddelta*(2.*pow(a, 2) - r*(2. + la*r + 2.*s));

	return dUR;
}

double potential_dUI(double r, hbl_parameters params){
	double a = params.a, om = params.om;
	int s = params.s, H = params.H, m = params.m;
	double delta = r*r - 2.*r + a*a;
	double ddelta = 2.*(r - 1.);
	double dUI;

	dUI = -2.*a*delta*(m + a*H*om) - 2.*a*ddelta*(m + a*H*om)*r
		+ 2.*om*pow(r, 2)*(delta*(1. - H)
		+ ddelta*(1. - H)*r - 2.*(1. + H)*r)*s + 4.*om*r*(delta*(1. - H)*r
		- (1. + H)*(-pow(a, 2) + pow(r, 2)))*s;

	return dUI;
}

double rstar_to_r(double rstar, double a){
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type* T;
	gsl_root_fsolver* s;
	double r = 0;
	double rp = 1 + sqrt(1-a*a);

	double x_lo = rp + 1.0e-10, x_hi = rstar;
	if(rstar < 5.0){
		x_hi = 5.0;
	}

	gsl_function F;
	rstar_parameters params = {rstar, a};

	F.function = &r_to_rstar_root;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do{
		iter++;
		status = gsl_root_fsolver_iterate(s);
	} while(status == GSL_CONTINUE && iter < max_iter);
	r = gsl_root_fsolver_root(s);

	gsl_root_fsolver_free(s);

	return r;
}

double r_to_rstar_root(double r, void* parameters){
	rstar_parameters params = *(rstar_parameters *)parameters;
	double a = params.a, rstar = params.rstar;

	return r_to_rstar(r, a) - rstar;
}



/////////////////////////////////////////
// Asymptotic series at the boundaries //
/////////////////////////////////////////

Result teuk_up_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r){
	double lambda = swsh_eigenvalue(s, L, m, a*omega);
	return teuk_up_asymptotic_infinity(a, s, L, m, omega, lambda, r);
}

Result teuk_up_asymptotic_infinity(const double &a, const int &s, const int &, const int &m, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - a*a);
	double rm = 1. - kappa, rp = 1. + kappa;
	double tau = (2.*omega - m*a)/kappa;
	double epsilonPlus = omega + 0.5*tau, epsilonMinus = omega - 0.5*tau;
	Complex cs = Complex(s);

	// aCH and bCH are the parameters that define the transformation between the Teukolsky equation and
	// the confluent Heun equation
	Complex aCH = I*epsilonMinus; // valid choices are aCH = - s - I*epsilonMinus or I*epsilonMinus
	Complex bCH = I*epsilonPlus; // valid choices are bCH = - s - I*epsilonPlus or I*epsilonPlus

	// the confluent Heun parameters
	Complex gammaCH = 1. + cs + 2.*aCH;
	Complex deltaCH = 1. + cs + 2.*bCH;
	Complex epsilonCH = 4.*I*omega*kappa;
	Complex alphaCH = epsilonCH*(1. + 2.*cs - 2.*I*omega + aCH + bCH);
	Complex qCH = -(bCH + aCH)*(1. + cs) - 2.*aCH*bCH + lambda - 2.*epsilonPlus*epsilonMinus + 2.*omega*Complex(m)*a
		- 2.*omega*(2*omega + I*cs)*(1. - kappa) + 2.*I*omega*cs*kappa + 2.*I*omega*kappa*(1. + 2.*aCH);

	// Asymptotic amplitude chosen to match normalization of MST solutions
	Complex Ctrans = pow(2., -2.*I*omega);
	Complex prefactor = pow((r - rp)/(r - rm), bCH)*pow(r - rm, -2.*cs - 1. + 2.*I*omega)*exp(I*omega*r);

	// recurrence relation for the asymptotic series coefficients for the confluent Heun solution near infinity
	int k = 0;
	Complex term = 1.;
	Complex sum = term;
	Complex maxTerm = term;
	k++;

	Complex cm2 = 0.;
	Complex cm1 = 1.;
	Complex c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
	c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
		Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
	c0 /= Complex(k)*pow(epsilonCH, 3);
	Complex arg = (2.*kappa)/(r - rm);

	term = c0*pow(arg, k);
	sum += term;
	maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
	k++;
	Complex previousTerm = term;
	while( (std::abs(term/sum) > DBL_EPSILON && std::abs(previousTerm/sum) > DBL_EPSILON && k < 100) || k <= 6 ){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);

		previousTerm = term;
		term = c0*pow(arg, k);
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		k++;
		// std::cout << "term = "<< term << ", sum = " << sum << "\n";
	}

	double error = std::abs(maxTerm/sum)*DBL_EPSILON;
	error = (std::abs(term/sum) < error)? error : std::abs(term/sum);

	return Result(Ctrans*prefactor*sum, Ctrans*prefactor*sum*error);
}

Result teuk_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r){
	double lambda = swsh_eigenvalue(s, L, m, a*omega);
	return teuk_up_derivative_asymptotic_infinity(a, s, L, m, omega, lambda, r);
}

Result teuk_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &, const int &m, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - a*a);
	double rm = 1. - kappa, rp = 1. + kappa;
	double tau = (2.*omega - m*a)/kappa;
	double epsilonPlus = omega + tau/2., epsilonMinus = omega - tau/2.;
	Complex cs = Complex(s);

	// aCH and bCH are the parameters that define the transformation between the Teukolsky equation and
	// the confluent Heun equation
	Complex aCH = I*epsilonMinus; // valid choices are aCH = - s - I*epsilonMinus or I*epsilonMinus
	Complex bCH = I*epsilonPlus; // valid choices are bCH = - s - I*epsilonPlus or I*epsilonPlus

	// the confluent Heun parameters
	Complex gammaCH = 1. + cs + 2.*aCH;
	Complex deltaCH = 1. + cs + 2.*bCH;
	Complex epsilonCH = 4.*I*omega*kappa;
	Complex alphaCH = epsilonCH*(1. + 2.*cs - 2.*I*omega + aCH + bCH);
	Complex qCH = -(bCH + aCH)*(1. + cs) - 2.*aCH*bCH + lambda - 2.*epsilonPlus*epsilonMinus + 2.*omega*Complex(m)*a
		- 2.*omega*(2*omega + I*cs)*(1. - kappa) + 2.*I*omega*cs*kappa + 2.*I*omega*kappa*(1. + 2.*aCH);

	// recurrence relation for the asymptotic series coefficients for the confluent Heun solution near infinity
	int k = 0;
	Complex gr = bCH/(r - rp) - (bCH + 2.*cs + 1. - 2.*I*omega)/(r - rm) + I*omega;
	Complex arg = (2.*kappa)/(r - rm);
	Complex term = gr;
	Complex sum = term;
	Complex maxTerm = term;
	k++;

	Complex cm2 = 0.;
	Complex cm1 = 1.;
	Complex c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
	c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
		Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
	c0 /= Complex(k)*pow(epsilonCH, 3);

	term = c0*(gr - Complex(k)/(r - rm))*pow(arg, k);
	sum += term;
	maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
	k++;
	Complex previousTerm = term;
	while( (std::abs(term/sum) > DBL_EPSILON && std::abs(previousTerm/sum) > DBL_EPSILON && k < 100 ) || k <= 5 ){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);

		previousTerm = term;
		term = c0*(gr - Complex(k)/(r - rm))*pow(arg, k);
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		k++;
	}

	double error = std::abs(maxTerm/sum)*DBL_EPSILON;
	error = (std::abs(term/sum) < error)? error : std::abs(term/sum);

	// Asymptotic amplitude chosen to match normalization of MST solutions
	Complex Ctrans = pow(2., -2.*I*omega);
	Complex prefactor = pow((r - rp)/(r - rm), bCH)*pow(r - rm, -2.*cs - 1. + 2.*I*omega)*exp(I*omega*r);

	return Result(Ctrans*prefactor*sum, Ctrans*prefactor*sum*error);
}

Result teuk_in_asymptotic_horizon(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r){
	double lambda = swsh_eigenvalue(s, L, m, a*omega);
	return teuk_in_asymptotic_horizon(a, s, L, m, omega, lambda, r);
}

Result teuk_in_asymptotic_horizon(const double &a, const int &s, const int &, const int &m, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - a*a);
	double rm = 1. - kappa, rp = 1. + kappa;
	double tau = (2.*omega - m*a)/kappa;
	double epsilonPlus = omega + tau/2., epsilonMinus = omega - tau/2.;
	double kfreq = omega - m*a/(2.*rp);
	Complex cs = Complex(s);

	// aCH and bCH are the parameters that define the transformation between the Teukolsky equation and
	// the confluent Heun equation
	double sgn = -1.; // another valid choice is sgn = +1.;
	Complex aCH = -cs - I*epsilonMinus; // another valid choice is aCH = I*epsilonMinus
	Complex bCH = I*epsilonPlus;

	// the confluent Heun parameters
	Complex gammaCH = 1. + cs + 2.*aCH;
	Complex deltaCH = 1. + cs + 2.*bCH;
	Complex epsilonCH = sgn*4.*I*omega*kappa;
	Complex alphaCH = epsilonCH*(1. + sgn*(cs - 2.*I*omega) + cs + aCH + bCH);
	Complex qCH = -(bCH + aCH)*(1. + cs) - 2.*aCH*bCH + lambda - 2.*epsilonPlus*epsilonMinus + 2.*omega*Complex(m)*a
		- 2.*omega*(2*omega + I*cs)*(1. - kappa) + sgn*2.*I*omega*cs*kappa + sgn*2.*I*omega*kappa*(1. + 2.*aCH);

	// recurrence relation for the asymptotic series coefficients for the confluent Heun solution near infinity
	int k = 0;
	Complex term = 1.;
	Complex sum = term;
	Complex maxTerm = term;
	k++;

	Complex cm2 = 0.;
	Complex cm1 = 1.;
	Complex c0 = -(alphaCH + epsilonCH*(Complex(k) - deltaCH - 1.))*cm2;
	c0 += -(Complex(k)*Complex(k) + Complex(k)*(gammaCH - deltaCH + epsilonCH - 1.) - deltaCH*(gammaCH + epsilonCH - 1.)
		+ alphaCH - qCH)*cm1;
	c0 /= Complex(k)*(Complex(k) - deltaCH + 1.);
	Complex arg = (r - rp)/(2.*kappa);

	term = c0*pow(arg, k);
	sum += term;
	maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
	k++;
	Complex previousTerm = term;
	while( (std::abs(term/sum) > DBL_EPSILON && std::abs(previousTerm/sum) > DBL_EPSILON && k < 300) || k <= 4 ){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - deltaCH - 1.))*cm2;
		c0 += -(Complex(k)*Complex(k) + Complex(k)*(gammaCH - deltaCH + epsilonCH - 1.) - deltaCH*(gammaCH + epsilonCH - 1.)
			+ alphaCH - qCH)*cm1;
		c0 /= Complex(k)*(Complex(k) - deltaCH + 1.);

		previousTerm = term;
		term = c0*pow(arg, k);
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		//std::cout << "sum_" << k << " = " << sum << " \n";
		k++;
	}

	double error = std::abs(maxTerm/sum)*DBL_EPSILON;
	error = (std::abs(term/sum) < error)? error : std::abs(term/sum);

	// Asymptotic amplitude chosen to match normalization of MST solutions
	Complex Btrans = pow(2., I*kfreq*rp/kappa - cs)*pow(kappa, I*kfreq*rm/kappa - cs)*exp(-I*kfreq*rp);
	Complex prefactor = pow((r - rm)/(2.*kappa), aCH)*pow(r - rp, -cs - bCH)*exp(sgn*I*omega*(r - rp));

	return Result(Btrans*prefactor*sum, Btrans*prefactor*sum*error);
}

Result teuk_in_derivative_asymptotic_horizon(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r){
	double lambda = swsh_eigenvalue(s, L, m, a*omega);
	return teuk_in_derivative_asymptotic_horizon(a, s, L, m, omega, lambda, r);
}

Result teuk_in_derivative_asymptotic_horizon(const double &a, const int &s, const int &, const int &m, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - a*a);
	double rm = 1. - kappa, rp = 1. + kappa;
	double tau = (2.*omega - m*a)/kappa;
	double epsilonPlus = omega + tau/2., epsilonMinus = omega - tau/2.;
	double kfreq = omega - m*a/(2.*rp);
	Complex cs = Complex(s);

	// aCH and bCH are the parameters that define the transformation between the Teukolsky equation and
	// the confluent Heun equation
	double sgn = -1.; // another valid choice is sgn = +1.;
	Complex aCH = -cs - I*epsilonMinus; // another valid choice is aCH = I*epsilonMinus
	Complex bCH = I*epsilonPlus;

	// the confluent Heun parameters
	Complex gammaCH = 1. + cs + 2.*aCH;
	Complex deltaCH = 1. + cs + 2.*bCH;
	Complex epsilonCH = sgn*4.*I*omega*kappa;
	Complex alphaCH = epsilonCH*(1. + sgn*(cs - 2.*I*omega) + cs + aCH + bCH);
	Complex qCH = -(bCH + aCH)*(1. + cs) - 2.*aCH*bCH + lambda - 2.*epsilonPlus*epsilonMinus + 2.*omega*Complex(m)*a
		- 2.*omega*(2*omega + I*cs)*(1. - kappa) + sgn*2.*I*omega*cs*kappa + sgn*2.*I*omega*kappa*(1. + 2.*aCH);

	// recurrence relation for the asymptotic series coefficients for the confluent Heun solution near infinity
	int k = 0;
	Complex gr = aCH/(r - rm) - (bCH + cs)/(r - rp) + sgn*I*omega;
	Complex arg = (r - rp)/(2.*kappa);
	Complex term = gr;
	Complex sum = term;
	Complex maxTerm = term;
	k++;

	Complex cm2 = 0.;
	Complex cm1 = 1.;
	Complex c0 = -(alphaCH + epsilonCH*(Complex(k) - deltaCH - 1.))*cm2;
	c0 += -(Complex(k)*Complex(k) + Complex(k)*(gammaCH - deltaCH + epsilonCH - 1.) - deltaCH*(gammaCH + epsilonCH - 1.)
		+ alphaCH - qCH)*cm1;
	c0 /= Complex(k)*(Complex(k) - deltaCH + 1.);

	term = c0*(gr*pow(arg, k) + Complex(k)/(2.*kappa)*pow(arg, k-1));
	sum += term;
	maxTerm = term;
	k++;

	Complex previousTerm = term;
	while( (std::abs(term/sum) > DBL_EPSILON && std::abs(previousTerm/sum) > DBL_EPSILON && k < 300) || k <= 4 ){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - deltaCH - 1.))*cm2;
		c0 += -(Complex(k)*Complex(k) + Complex(k)*(gammaCH - deltaCH + epsilonCH - 1.) - deltaCH*(gammaCH + epsilonCH - 1.)
			+ alphaCH - qCH)*cm1;
		c0 /= Complex(k)*(Complex(k) - deltaCH + 1.);

		previousTerm = term;
		term = c0*(gr*pow(arg, k) + Complex(k)/(2.*kappa)*pow(arg, k-1));
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		k++;
	}

	double error = std::abs(maxTerm/sum)*DBL_EPSILON;
	error = (std::abs(term/sum) < error)? error : std::abs(term/sum);

	// Asymptotic amplitude chosen to match normalization of MST solutions
	Complex Btrans = pow(2., I*kfreq*rp/kappa - cs)*pow(kappa, I*kfreq*rm/kappa - cs)*exp(-I*kfreq*rp);
	Complex prefactor = pow((r - rm)/(2.*kappa), aCH)*pow(r - rp, -cs - bCH)*exp(sgn*I*omega*(r - rp));

	return Result(Btrans*prefactor*sum, Btrans*prefactor*sum*error);
}

Result gsn_up_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r){
	double lambda = swsh_eigenvalue(s, L, m, a*omega);
	return gsn_up_asymptotic_infinity(a, s, L, m, omega, lambda, r);
}

Result gsn_up_asymptotic_infinity(const double &r, hbl_parameters params){
	return gsn_up_asymptotic_infinity(params.a, params.s, params.m, params.om, params.la, r);
}

Complex gsn_asymptotic_initial_sum(const double &r, hbl_parameters params){
	return gsn_asymptotic_initial_sum(params.a, params.s, params.m, params.om, params.la, r);
}

Complex gsn_asymptotic_initial_sum(const double &a, const int &, const int &m, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - pow(a, 2));
	if(std::abs(a) < DBL_EPSILON){
		return (Complex(0.,-0.02083333333333333333)*pow(omega,-3)*pow(r,-2)*(-8. + 12.*r -
			6.*pow(r,2) + pow(r,3))*(r*(2.*lambda - Complex(0.,12.)*omega + pow(lambda,2))*(-12.*(2. + lambda)
				+ r*(12. + 8.*lambda + Complex(0.,36.)*omega + pow(lambda,2)) - Complex(0.,12.)*omega*pow(r,2)))*pow((-2. + r)*r,-3));
	}
	return (Complex(0.,-0.020833333333333333333)*pow(omega,-3)*pow(r,-2)*(-4.*(1. + kappa) + 6.*(1. + kappa)*r + (3. + kappa - 3.*r)*pow(a,2) -
		3.*(1. + kappa)*pow(r,2) + pow(r,3))*(r*(2.*lambda - Complex(0.,12.)*omega + pow(lambda,2))*(-12.*(2. + lambda)
			+ r*(12. + 8.*lambda + Complex(0.,36.)*omega + pow(lambda,2)) - Complex(0.,12.)*omega*pow(r,2)) + 8.*a*m*omega*r*(-6.*(6. + 7.*lambda) +
				r*(18. + 27.*lambda + Complex(0.,36.)*omega + 5.*pow(lambda,2)) - Complex(0.,18.)*omega*pow(r,2)) +
					24.*pow(a,4)*pow(omega,2)*(6. - 7.*lambda + Complex(0.,6.)*omega*r + 6.*pow(omega,2)*pow(r,2)) -
						24.*m*omega*pow(a,3)*(-6. - 7.*lambda + Complex(0.,6.)*omega*r + 12.*pow(omega,2)*pow(r,2)) +
							2.*pow(a,2)*(3.*pow(lambda,3) - 12.*lambda*(-1. + Complex(0.,1.)*omega*(3. + r) + r*(-14. + 5.*r)*pow(omega,2)) -
								2.*pow(lambda,2)*(-6. + Complex(0.,3.)*omega*r + 10.*pow(omega,2)*pow(r,2)) + 72.*omega*(Complex(0.,-1.) - 3.*omega*r +
									omega*(1. - Complex(0.,2.)*omega + pow(m,2))*pow(r,2) + Complex(0.,1.)*pow(omega,2)*pow(r,3))))*pow((-2. + r)*r + pow(a,2),-3));
}

Result gsn_up_asymptotic_infinity(const double &a, const int &s, const int &, const int &m, const double &omega, const double &lambda, const double &r){
	return gsn_up_asymptotic_infinity(a, s, m, omega, lambda, r);
}

Result gsn_up_asymptotic_infinity(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r){
	Result chi = gsn_up_asymptotic_infinity_chi_series(a, s, m, omega, lambda, r);
	hbl_parameters params = {
		.a = a,
		.s = s,
		.m = m,
		.om = omega,
		.la = lambda,
		.H = 0
	};
	return GSN_chi_to_GSN_X(chi, r, params);
}

Result gsn_up_asymptotic_infinity_chi_series(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - a*a);
	double rm = 1. - kappa, rp = 1. + kappa;
	double tau = (2.*omega - m*a)/kappa;
	double epsilonPlus = omega + tau/2., epsilonMinus = omega - tau/2.;
	if(double(a) < DBL_EPSILON){
		kappa = 1.;
		rm = 0.;
		rp = 2.;
		tau = 2.*omega;
		epsilonPlus = 2.*omega;
		epsilonMinus = 0.;
	}
	Complex cs = Complex(s);
	hbl_parameters params = {
		.a = a,
		.s = s,
		.m = m,
		.om = omega,
		.la = lambda,
		.H = 0
	};

	// aCH and bCH are the parameters that define the transformation between the Teukolsky equation and
	// the confluent Heun equation
	Complex aCH = -cs - I*epsilonMinus;
	Complex bCH = I*epsilonPlus;

	// Asymptotic amplitude chosen to match normalization of MST solutions
	Complex Ctrans = pow(2., -2.*I*omega);
	Complex prefactor = pow((r - rp)/(r - rm), bCH)*pow(r - rm, -2.*cs - 1. + 2.*I*omega)*exp(I*omega*r);

	// the confluent Heun parameters
	Complex gammaCH = 1. + cs + 2.*aCH;
	Complex deltaCH = 1. + cs + 2.*bCH;
	Complex epsilonCH = 4.*I*omega*kappa;
	Complex alphaCH = epsilonCH*(1. + 2.*cs - 2.*I*omega + aCH + bCH);
	Complex qCH = -(bCH + aCH)*(1. + cs) - 2.*aCH*bCH + lambda - 2.*epsilonPlus*epsilonMinus + 2.*omega*Complex(m)*a
		- 2.*omega*(2*omega + I*cs)*(1. - kappa) + 2.*I*omega*cs*kappa + 2.*I*omega*kappa*(1. + 2.*aCH);
	if(double(a) < DBL_EPSILON){
		qCH = -(bCH + aCH)*(1. + cs) - 2.*aCH*bCH + lambda + 2.*I*omega*(1. + 2.*aCH + cs);
	}

	// recurrence relation for the asymptotic series coefficients for the confluent Heun solution near infinity
	int k = 0;
	Complex alphaGSN = GSN_alpha(r, params), betaGSN = GSN_beta(r, params)*pow(teuk_delta(r, params), s + 1);
	Complex gr = bCH/(r - rp) - (bCH + 2.*cs + 1. - 2.*I*omega)/(r - rm) + I*omega;
	Complex arg = (2.*kappa)/(r - rm);
	Complex alphaBetaGSN = (alphaGSN + betaGSN*gr);

	Complex term = gsn_asymptotic_initial_sum(r, params);
	Complex sum = term;
	Complex maxTerm = term;
	k++;

	Complex cm2 = 0.;
	Complex cm1 = 1.;
	Complex c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
	c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
		Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
	c0 /= Complex(k)*pow(epsilonCH, 3);

	for(k = 2; k < 5; k++){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);
	}

	k = 4;
	term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
	sum += term;
	maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
	k++;
	Complex previousTerm = term;
	while( (std::abs(term/sum) > DBL_EPSILON || std::abs(previousTerm/sum) > DBL_EPSILON) && k < 300){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);

		previousTerm = term;
		term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		k++;
	}
	cm2 = cm1;
	cm1 = c0;
	c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
	c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
		Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
	c0 /= Complex(k)*pow(epsilonCH, 3);

	previousTerm = term;
	term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
	sum += term;
	maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
	k++;
	while( (std::abs(term/sum) > DBL_EPSILON || std::abs(previousTerm/sum) > DBL_EPSILON) && k < 300){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);

		previousTerm = term;
		term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		k++;
	}

	double error = std::abs(maxTerm/sum)*DBL_EPSILON;
	error = (std::abs(term/sum) < error)? error : std::abs(term/sum);

	return Result(Ctrans*prefactor*sum, Ctrans*prefactor*sum*error);
}

Result gsn_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &L, const int &m, const double &omega, const double &r){
	double lambda = swsh_eigenvalue(s, L, m, a*omega);
	return gsn_up_derivative_asymptotic_infinity(a, s, m, omega, lambda, r);
}

Result gsn_up_derivative_asymptotic_infinity(const double &r, hbl_parameters params){
	return gsn_up_derivative_asymptotic_infinity(params.a, params.s, params.m, params.om, params.la, r);
}

Complex gsn_asymptotic_derivative_initial_sum_series(const double &r, hbl_parameters params){
	return gsn_asymptotic_derivative_initial_sum_series(params.a, params.s, params.m, params.om, params.la, r);
}

Complex gsn_asymptotic_derivative_initial_sum_series(const double &a, const int &s, const int &mTemp, const double &omega, const double &lambda, const double &r){
	Complex sum = 0.;
	sum += gsn_asymptotic_derivative_initial_sum_term_1(a, s, mTemp, omega, lambda, r)*pow(omega*r, -1);
	sum += gsn_asymptotic_derivative_initial_sum_term_2(a, s, mTemp, omega, lambda, r)*pow(omega*r, -2);
	sum += gsn_asymptotic_derivative_initial_sum_term_3(a, s, mTemp, omega, lambda, r)*pow(omega*r, -3);
	sum += gsn_asymptotic_derivative_initial_sum_term_4(a, s, mTemp, omega, lambda, r)*pow(omega*r, -4);
	sum += gsn_asymptotic_derivative_initial_sum_term_5(a, s, mTemp, omega, lambda, r)*pow(omega*r, -5);
	sum += gsn_asymptotic_derivative_initial_sum_term_6(a, s, mTemp, omega, lambda, r)*pow(omega*r, -6);
	sum += gsn_asymptotic_derivative_initial_sum_term_7(a, s, mTemp, omega, lambda, r)*pow(omega*r, -7);
	sum += gsn_asymptotic_derivative_initial_sum_term_8(a, s, mTemp, omega, lambda, r)*pow(omega*r, -8);
	sum += gsn_asymptotic_derivative_initial_sum_term_9(a, s, mTemp, omega, lambda, r)*pow(omega*r, -9);
	sum += gsn_asymptotic_derivative_initial_sum_term_10(a, s, mTemp, omega, lambda, r)*pow(omega*r, -10);
	sum += gsn_asymptotic_derivative_initial_sum_term_11(a, s, mTemp, omega, lambda, r)*pow(omega*r, -11);
	sum += gsn_asymptotic_derivative_initial_sum_term_12(a, s, mTemp, omega, lambda, r)*pow(omega*r, -12);
	return sum;
}

Complex gsn_asymptotic_derivative_initial_sum(const double &r, hbl_parameters params){
	return gsn_asymptotic_derivative_initial_sum(params.a, params.s, params.m, params.om, params.la, r);
}

Complex gsn_asymptotic_derivative_initial_sum(const double &a, const int &, const int &mTemp, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - pow(a, 2));
	if(std::abs(a) < DBL_EPSILON){
		kappa = 1.;
	}
	Complex m = Complex(mTemp);
	return 0.0026041666666666665*pow(omega, -4)*pow(r,-3)*pow((-2. + r)*r + pow(a,2),-6)*((48.*a*lambda*omega*(m + a*omega) + 4.*pow(lambda,3) + pow(lambda,4) + pow(lambda,2)*(4. + 40.*a*m*omega - 40.*pow(a,2)*pow(omega,2)) + 144.*pow(omega,2)*(1. - 2.*m*omega*pow(a,3) + pow(a,2)*pow(m,2) + pow(a,4)*pow(omega,2)))*pow(1. + kappa - 1.*r,4)*(Complex(0.,-6.)*m*r*pow(a,5) + Complex(0.,6.)*(Complex(0.,2.) + omega*r)*pow(a,6) + m*pow(a,3)*pow(r,2)*(Complex(0.,-4.)*(-5. + kappa) + Complex(0.,1.)*lambda*r + 4.*(Complex(0.,-2.) + omega + kappa*omega)*r - 10.*omega*pow(r,2)) + a*m*pow(r,3)*(Complex(0.,2.)*(-11. + kappa) + Complex(0.,2.)*(12. + kappa - 1.*lambda)*r + (Complex(0.,-8.) + Complex(0.,1.)*lambda + 4.*(4. + kappa)*omega)*pow(r,2) - 10.*omega*pow(r,3)) + r*pow(a,4)*(10.*(7. + kappa) - 2.*(17. + 2.*lambda - Complex(0.,2.)*(-5. + kappa)*omega)*r - 1.*omega*(Complex(0.,1.)*lambda + 2.*(Complex(0.,2.) + omega + kappa*omega))*pow(r,2) + 8.*pow(omega,2)*pow(r,3)) + pow(a,2)*pow(r,2)*(-24.*(5. + kappa) + r*(102. + 15.*lambda + kappa*(6. + lambda - Complex(0.,2.)*omega) + Complex(0.,22.)*omega - 2.*(1. + kappa)*pow(m,2)) + (-18. - 7.*lambda + Complex(0.,18.)*omega + Complex(0.,4.)*kappa*omega + Complex(0.,2.)*lambda*omega + 2.*pow(m,2))*pow(r,2) - 2.*omega*(Complex(0.,5.) + Complex(0.,1.)*lambda + 2.*(4. + kappa)*omega)*pow(r,3) + 16.*pow(omega,2)*pow(r,4)) + pow(r,3)*(12.*(5. + kappa) - 2.*(32. + 7.*lambda + kappa*(2. + lambda))*r + (16. + (13. + kappa)*lambda - Complex(0.,6.)*(5. + kappa)*omega)*pow(r,2) + Complex(0.,1.)*(18.*omega + lambda*(Complex(0.,3.) + 2.*omega))*pow(r,3) - 1.*omega*(Complex(0.,1.)*lambda + 2.*(7. + kappa)*omega)*pow(r,4) + 8.*pow(omega,2)*pow(r,5))) - 8.*omega*((Complex(0.,4.) - 6.*(-1. + kappa)*omega)*pow(lambda,2) + Complex(0.,1.)*pow(lambda,3) + 4.*lambda*(Complex(0.,1.) + (6. - 3.*kappa + Complex(0.,7.)*a*m)*omega - Complex(0.,1.)*(12. - 12.*kappa + pow(a,2))*pow(omega,2)) + 24.*omega*(1. + a*m*(Complex(0.,1.) - 3.*(-1. + kappa)*omega) + Complex(0.,1.)*omega*(-3. + 3.*kappa + pow(a,2)) + (8.*(-1. + kappa) + (3. + kappa)*pow(a,2))*pow(omega,2)))*pow(1. + kappa - 1.*r,3)*(Complex(0.,6.)*m*r*pow(a,3) + 6.*(2. - Complex(0.,1.)*omega*r)*pow(a,4) + a*m*pow(r,2)*(Complex(0.,-12.) - Complex(0.,1.)*(-6. + lambda)*r + 6.*omega*pow(r,2)) + r*pow(a,2)*(-36. + 4.*(3. + lambda + Complex(0.,3.)*omega)*r + Complex(0.,1.)*(6. + lambda)*omega*pow(r,2) - 6.*pow(omega,2)*pow(r,3)) + pow(r,2)*(24. - 6.*(2. + lambda)*r + 2.*(lambda - Complex(0.,6.)*omega)*pow(r,2) + Complex(0.,1.)*lambda*omega*pow(r,3) - 6.*pow(omega,2)*pow(r,4)))*pow((-2. + r)*r + pow(a,2),2) + 48.*(lambda*(2. + Complex(0.,8.)*(-1. + kappa)*omega) + 12.*omega*(Complex(0.,-1.) + a*m + omega*(-4. + 4.*kappa + pow(a,2))) + pow(lambda,2))*pow(omega,2)*pow(1. + kappa - 1.*r,2)*(Complex(0.,6.)*m*r*pow(a,5) + 6.*(2. - Complex(0.,1.)*omega*r)*pow(a,6) + m*pow(a,3)*pow(r,2)*(Complex(0.,-4.)*(7. + kappa) - Complex(0.,1.)*lambda*r + 4.*(Complex(0.,4.) + omega + kappa*omega)*r + 2.*omega*pow(r,2)) + a*m*pow(r,3)*(Complex(0.,2.)*(13. + kappa) + Complex(0.,2.)*(-12. + kappa + lambda)*r + (Complex(0.,4.) - Complex(0.,1.)*lambda + 4.*(-2. + kappa)*omega)*pow(r,2) + 2.*omega*pow(r,3)) + r*pow(a,4)*(10.*(-5. + kappa) + 2.*(7. + 2.*lambda + Complex(0.,2.)*(7. + kappa)*omega)*r - 1.*omega*(Complex(0.,-1.)*lambda + 2.*(Complex(0.,2.) + omega + kappa*omega))*pow(r,2) - 4.*pow(omega,2)*pow(r,3)) + pow(a,2)*pow(r,2)*(-24.*(-3. + kappa) - 1.*r*(42. + 13.*lambda - 1.*kappa*(6. + lambda - Complex(0.,2.)*omega) + Complex(0.,26.)*omega + 2.*(1. + kappa)*pow(m,2)) + (6. + lambda*(5. - Complex(0.,2.)*omega) - Complex(0.,6.)*omega + Complex(0.,4.)*kappa*omega + 2.*pow(m,2))*pow(r,2) + 2.*omega*(Complex(0.,1.) + Complex(0.,1.)*lambda - 2.*(-2. + kappa)*omega)*pow(r,3) - 8.*pow(omega,2)*pow(r,4)) + pow(r,3)*(12.*(-3. + kappa) - 2.*(-16. - 5.*lambda + kappa*(2. + lambda))*r + (-8. + (-7. + kappa)*lambda - Complex(0.,6.)*(-3. + kappa)*omega)*pow(r,2) + (lambda - Complex(0.,6.)*omega - Complex(0.,2.)*lambda*omega)*pow(r,3) + omega*(Complex(0.,1.)*lambda - 2.*(-5. + kappa)*omega)*pow(r,4) - 4.*pow(omega,2)*pow(r,5)))*pow((-2. + r)*r + pow(a,2),2) + 192.*(Complex(0.,1.)*lambda - 6.*(-1. + kappa)*omega)*(1. + kappa - 1.*r)*pow(omega,3)*(Complex(0.,6.)*m*r*pow(a,5) + 6.*(2. - Complex(0.,1.)*omega*r)*pow(a,6) + m*pow(a,3)*pow(r,2)*(Complex(0.,-8.)*(4. + kappa) + (Complex(0.,20.) - Complex(0.,1.)*lambda + 8.*(1. + kappa)*omega)*r - 2.*omega*pow(r,2)) + a*m*pow(r,3)*(Complex(0.,4.)*(7. + kappa) + Complex(0.,2.)*(-12. + 2.*kappa + lambda)*r + (Complex(0.,2.) - Complex(0.,1.)*lambda - 4.*omega + 8.*kappa*omega)*pow(r,2) - 2.*omega*pow(r,3)) + r*pow(a,4)*(20.*(-2. + kappa) + 4.*(1. + lambda + Complex(0.,2.)*(4. + kappa)*omega)*r - 1.*omega*(Complex(0.,-1.)*lambda + 4.*(Complex(0.,2.) + omega + kappa*omega))*pow(r,2) - 2.*pow(omega,2)*pow(r,3)) - 2.*pow(a,2)*pow(r,2)*(24.*(-1. + kappa) + r*(6. + 6.*lambda - 1.*kappa*(6. + lambda - Complex(0.,2.)*omega) + Complex(0.,14.)*omega + 2.*(1. + kappa)*pow(m,2)) - 1.*(lambda*(2. - Complex(0.,1.)*omega) + Complex(0.,4.)*kappa*omega + 2.*pow(m,2))*pow(r,2) + omega*(Complex(0.,1.) - Complex(0.,1.)*lambda + (-2. + 4.*kappa)*omega)*pow(r,3) + 2.*pow(omega,2)*pow(r,4)) - 1.*pow(r,3)*(-24.*(-1. + kappa) + 4.*(-2. + kappa)*(2. + lambda)*r - 2.*(-2. + (-2. + kappa)*lambda - Complex(0.,6.)*(-1. + kappa)*omega)*pow(r,2) + Complex(0.,2.)*lambda*omega*pow(r,3) + omega*(Complex(0.,-1.)*lambda + 4.*(-2. + kappa)*omega)*pow(r,4) + 2.*pow(omega,2)*pow(r,5)))*pow((-2. + r)*r + pow(a,2),3) + 384.*pow(omega,4)*(Complex(0.,-6.)*m*r*pow(a,5) + Complex(0.,6.)*(Complex(0.,2.) + omega*r)*pow(a,6) + m*pow(a,3)*pow(r,2)*(Complex(0.,12.)*(3. + kappa) + Complex(0.,1.)*lambda*r - 12.*(Complex(0.,2.) + omega + kappa*omega)*r + 6.*omega*pow(r,2)) + r*pow(a,4)*(-30.*(-1. + kappa) - 2.*(-3. + 2.*lambda + Complex(0.,6.)*(3. + kappa)*omega)*r + omega*(Complex(0.,-1.)*lambda + 6.*(Complex(0.,2.) + omega + kappa*omega))*pow(r,2)) + a*m*pow(r,3)*(Complex(0.,-6.)*(5. + kappa) - Complex(0.,2.)*(-12. + 3.*kappa + lambda)*r + Complex(0.,1.)*(lambda + Complex(0.,12.)*kappa*omega)*pow(r,2) + 6.*omega*pow(r,3)) + pow(a,2)*pow(r,2)*(-24. + 72.*kappa + r*(-18. + 11.*lambda - 3.*kappa*(6. + lambda - Complex(0.,2.)*omega) + Complex(0.,30.)*omega + 6.*(1. + kappa)*pow(m,2)) - 1.*(lambda*(3. - Complex(0.,2.)*omega) + Complex(0.,6.)*(Complex(0.,1.) + omega + 2.*kappa*omega) + 6.*pow(m,2))*pow(r,2) + 2.*omega*(Complex(0.,3.) - Complex(0.,1.)*lambda + 6.*kappa*omega)*pow(r,3)) + pow(r,3)*(12. - 6.*lambda*r + (lambda - Complex(0.,6.)*omega)*pow(r,2) + (lambda - Complex(0.,6.)*omega + Complex(0.,2.)*lambda*omega)*pow(r,3) + (Complex(0.,-1.)*lambda - 6.*omega)*omega*pow(r,4) + 3.*kappa*(-12. + 2.*(2. + lambda)*r - 1.*(lambda - Complex(0.,6.)*omega)*pow(r,2) + 2.*pow(omega,2)*pow(r,4))))*pow((-2. + r)*r + pow(a,2),4));
}

Result gsn_up_derivative_asymptotic_infinity(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r){
	Result chi = gsn_up_asymptotic_infinity_chi_series(a, s, m, omega, lambda, r);
	Result dchi = gsn_up_derivative_asymptotic_infinity_chi_series(a, s, m, omega, lambda, r);
	hbl_parameters params = {
		.a = a,
		.s = s,
		.m = m,
		.om = omega,
		.la = lambda,
		.H = 0
	};
	return GSN_dchi_to_GSN_dX(chi, dchi, r, params);
}

Result gsn_up_derivative_asymptotic_infinity_chi_series(const double &a, const int &s, const int &m, const double &omega, const double &lambda, const double &r){
	double kappa = sqrt(1. - a*a);
	double rm = 1. - kappa, rp = 1. + kappa;
	double tau = (2.*omega - m*a)/kappa;
	double epsilonPlus = omega + tau/2., epsilonMinus = omega - tau/2.;
	Complex cs = Complex(s);
	hbl_parameters params = {
		.a = a,
		.s = s,
		.m = m,
		.om = omega,
		.la = lambda,
		.H = 0
	};

	// aCH and bCH are the parameters that define the transformation between the Teukolsky equation and
	// the confluent Heun equation
	Complex aCH = -cs - I*epsilonMinus;
	Complex bCH = I*epsilonPlus;

	// Asymptotic amplitude chosen to match normalization of MST solutions
	Complex Ctrans = pow(2., -2.*I*omega);
	Complex prefactor = pow((r - rp)/(r - rm), bCH)*pow(r - rm, -2.*cs - 1. + 2.*I*omega)*exp(I*omega*r);
	// Complex prefactorP = I*omega + bCH*(rp - rm)/(r - rm)/(r - rp) + (2.*I*omega - 2.*s - 1.)/(r - rm);

	// the confluent Heun parameters
	Complex gammaCH = 1. + cs + 2.*aCH;
	Complex deltaCH = 1. + cs + 2.*bCH;
	Complex epsilonCH = 4.*I*omega*kappa;
	Complex alphaCH = epsilonCH*(1. + 2.*cs - 2.*I*omega + aCH + bCH);
	Complex qCH = -(bCH + aCH)*(1. + cs) - 2.*aCH*bCH + lambda - 2.*epsilonPlus*epsilonMinus + 2.*omega*Complex(m)*a
		- 2.*omega*(2*omega + I*cs)*(1. - kappa) + 2.*I*omega*cs*kappa + 2.*I*omega*kappa*(1. + 2.*aCH);

	// recurrence relation for the asymptotic series coefficients for the confluent Heun solution near infinity
	int k = 0;
	Complex alphaGSN = GSN_alphaP(r, params) + GSN_beta(r, params)*pow(teuk_delta(r, params), s)*GSN_V(r, params), betaGSN = GSN_alpha(r, params) + GSN_betaP(r, params)*pow(teuk_delta(r, params), s + 1);
	Complex gr = bCH/(r - rp) - (bCH + 2.*cs + 1. - 2.*I*omega)/(r - rm) + I*omega;
	Complex arg = (2.*kappa)/(r - rm);
	Complex alphaBetaGSN = (alphaGSN + betaGSN*gr);

	Complex term = gsn_asymptotic_derivative_initial_sum_series(r, params);
	Complex sum = term;
	Complex maxTerm = term;
	k++;

	Complex cm2 = 0.;
	Complex cm1 = 1.;
	Complex c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
	c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
		Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
	c0 /= Complex(k)*pow(epsilonCH, 3);

	for(k = 2; k < 6; k++){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);
	}

	k = 5;
	term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
	sum += term;
	maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
	k++;
	Complex previousTerm = term;
	while( (std::abs(term/sum) > DBL_EPSILON || std::abs(previousTerm/sum) > DBL_EPSILON) && k < 300){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);

		previousTerm = term;
		term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		k++;
	}
	cm2 = cm1;
	cm1 = c0;
	c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
	c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
		Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
	c0 /= Complex(k)*pow(epsilonCH, 3);

	previousTerm = term;
	term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
	sum += term;
	maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
	k++;
	while( (std::abs(term/sum) > DBL_EPSILON || std::abs(previousTerm/sum) > DBL_EPSILON) && k < 300){
		cm2 = cm1;
		cm1 = c0;
		c0 = -(alphaCH + epsilonCH*(Complex(k) - 2.))*(alphaCH + epsilonCH*(Complex(k) - gammaCH - 1.))*cm2;
		c0 += (alphaCH*alphaCH + alphaCH*epsilonCH*(2.*Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + epsilonCH*epsilonCH*(
			Complex(k)*(Complex(k) - gammaCH - deltaCH + epsilonCH - 1.) + gammaCH + deltaCH - epsilonCH - qCH))*cm1;
		c0 /= Complex(k)*pow(epsilonCH, 3);

		previousTerm = term;
		term = c0*(alphaBetaGSN - betaGSN*Complex(k)/(r - rm))*pow(arg, k);
		sum += term;
		maxTerm = (std::abs(term) < std::abs(maxTerm))? maxTerm : term;
		k++;
	}

	double error = std::abs(maxTerm/sum)*DBL_EPSILON;
	error = (std::abs(term/sum) < error)? error : std::abs(term/sum);

	Result term1 = Result(Ctrans*prefactor*sum, Ctrans*prefactor*sum*error);

	return term1;
}
