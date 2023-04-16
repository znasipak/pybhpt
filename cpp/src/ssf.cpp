// ssf.cpp

#include "ssf.hpp"

#define SSF_EPSILON 1.e-12
#define SSF_LMAX 40
#define SSF_KMAX 100
#define SSF_NMAX 200
#define ZERO_FREQ_LIMIT 1.e-9
#define CONVERGENCE_CRITERIA_LARGE 1000
#define COUPLING_EPSILON 1.e-24

double Clm(int l, int m){
  if(l >= abs(m) && l !=0 ){
    return sqrt(l + m)/sqrt(2.*l + 1.)*sqrt(l - m)/sqrt(2.*l - 1.);
  }else{
    return 0.;
  }
}

double delta_coupling(int coeff, int l, int m){
  if(coeff == 1){
    return l*Clm(l + 1, m);
  }else if(coeff == -1){
    return -(l + 1)*Clm(l, m);
  }else{
    return 0.;
  }
}

double zeta_coupling(int coeff, int l, int m){
  if(coeff == 1){
    return Clm(l + 1, m)*(l*(1. - pow(Clm(l + 1, m), 2) - pow(Clm(l + 2, m), 2)) + (l + 1.)*pow(Clm(l, m), 2));
  }else if(coeff == -1){
    return -Clm(l, m)*((l + 1.)*(1. - pow(Clm(l - 1, m), 2) - pow(Clm(l, m), 2)) + l*pow(Clm(l + 1, m), 2));
  }else if(coeff == -3){
    return (l + 1.)*Clm(l, m)*Clm(l - 1, m)*Clm(l - 2, m);
  }else if(coeff == 3){
    return -l*Clm(l + 1, m)*Clm(l + 2, m)*Clm(l + 3, m);
  }else{
    return 0.;
  }
}

// double beta_coupling(int coeff, int l, int m, double thp){
//   if(l >= abs(m)){
//     return 1.5*delta_coupling(coeff, l, m)/sin(thp) - 0.5*zeta_coupling(coeff, l, m)/(pow(sin(thp), 3));
//   }else{
//     return 0.;
//   }
// }

void beta_coupling(Vector &couplingVector, TeukolskyMode &teuk, int n, int l, int m){
	double bdelta = teuk.getCouplingCoefficient(l)*delta_coupling(n, l, m);
	double bzeta = teuk.getCouplingCoefficient(l)*zeta_coupling(n, l, m);
	for(size_t jth = 0; jth < couplingVector.size(); jth++){
		couplingVector[jth] += 1.5*bdelta/sin(teuk.getPolarPoints(jth)) - 0.5*bzeta/(pow(sin(teuk.getPolarPoints(jth)), 3));
	}
}

double beta_coupling(TeukolskyMode &teuk, int n, int l, int m, int jth){
	double bdelta = teuk.getCouplingCoefficient(l)*delta_coupling(n, l, m);
	double bzeta = teuk.getCouplingCoefficient(l)*zeta_coupling(n, l, m);
	return 1.5*bdelta/sin(teuk.getPolarPoints(jth)) - 0.5*bzeta/(pow(sin(teuk.getPolarPoints(jth)), 3));
}

double beta_coupling(int n, int l, int m, double theta){
	double bdelta = delta_coupling(n, l, m);
	double bzeta = zeta_coupling(n, l, m);
	return 1.5*bdelta/sin(theta) - 0.5*bzeta/(pow(sin(theta), 3));
}

// Coupling due to Warburton window function
void polar_coupling_coefficient_3(Vector &couplingVector, TeukolskyMode &teuk, int l, int m){
	for(int n = -3; n <= 3; n += 2){
		if(l - n >= abs(m)){
			beta_coupling(couplingVector, teuk, n, l - n, m);
		}
	}
}

double polar_coupling_coefficient_3(TeukolskyMode &teuk, int l, int m, int jth){
	double coupling = 0.;
	for(int n = -3; n <= 3; n += 2){
		if(l - n >= abs(m)){
			coupling += beta_coupling(teuk, n, l - n, m, jth);
		}
	}
	return coupling;
}

// Coupling due to Nasipak window function

double alpha_coupling(int n, int l, int m){
	switch(n){
		case -2:
			return -(l + 1)*Clm(l, m)*Clm(l - 1, m);
		case -1:
			return -(l + 1)*Clm(l, m);
		case 0:
			return l*pow(Clm(l + 1, m), 2) - (l + 1)*pow(Clm(l, m), 2);
		case 1:
			return l*Clm(l + 1, m);
		case 2:
			return l*Clm(l + 1, m)*Clm(l + 2, m);
		default:
			return 0.;
	}
}

double xi_coupling(int n, double thp){
	double cth = cos(thp);
	double sth = sin(thp);
  if(abs(n) == 2){
    return cth*pow(sth, -3);
  }else if(abs(n) == 1){
    return (1. - 2.*pow(cth, 2))*pow(sth, -3);
  }else if(abs(n) == 0){
    return cth*pow(sth, -3);
  }else{
    return 0.;
  }
}

void beta_coupling_nasipak(Vector &couplingVector, TeukolskyMode &teuk, int n, int l, int m){
	double balpha = teuk.getCouplingCoefficient(l)*alpha_coupling(n, l, m);
	for(size_t jth = 0; jth < couplingVector.size(); jth++){
		couplingVector[jth] += balpha*xi_coupling(n, teuk.getPolarPoints(jth));
	}
}

double beta_coupling_nasipak(TeukolskyMode &teuk, int n, int l, int m, int jth){
	double balpha = teuk.getCouplingCoefficient(l)*alpha_coupling(n, l, m);
	return balpha*xi_coupling(n, teuk.getPolarPoints(jth));
}

void polar_coupling_coefficient_2(Vector &couplingVector, TeukolskyMode &teuk, int l, int m){
	for(int n = -2; n <= 2; n++){
		if(l - n >= abs(m)){
			beta_coupling_nasipak(couplingVector, teuk, n, l - n, m);
		}
	}
}

double polar_coupling_coefficient_2(TeukolskyMode &teuk, int l, int m, int jth){
	double coupling = 0.;
	for(int n = -2; n <= 2; n++){
		if(l - n >= abs(m)){
			coupling += beta_coupling_nasipak(teuk, n, l - n, m, jth);
		}
	}
	return coupling;
}

double alpha_coupling_higher_order(int n, int k, int l, int m){
	if(n == -4 && k == 3){
		return -(l + 1)*Clm(l, m)*Clm(l - 1, m)*Clm(l - 2, m)*Clm(l - 3, m);
	}else if(n == -3 && k == 2){
		return -(l + 1)*Clm(l, m)*Clm(l - 1, m)*Clm(l - 2, m);
	}else if(n == -2 && k == 1){
		return -(l + 1)*Clm(l, m)*Clm(l - 1, m);
	}else if(n == -2 && k == 3){
		return Clm(l, m)*Clm(l - 1, m)*(l*pow(Clm(l + 1, m), 2) - (l + 1.)*(pow(Clm(l, m), 2) + pow(Clm(l - 1, m), 2) + pow(Clm(l - 2, m), 2)));
	}else if(n == -1 && k == 0){
		return -(l + 1)*Clm(l, m);
	}else if(n == -1 && k == 2){
		return Clm(l, m)*(l*pow(Clm(l + 1, m), 2) - (l + 1.)*(pow(Clm(l, m), 2) + pow(Clm(l - 1, m), 2)));
	}else if(n == 0 && k == 1){
		return l*pow(Clm(l + 1, m), 2) - (l + 1.)*pow(Clm(l, m), 2);
	}else if(n == 0 && k == 3){
		return l*pow(Clm(l + 1, m), 2)*(pow(Clm(l + 1, m), 2) + pow(Clm(l + 2, m), 2)) - (l + 1.)*pow(Clm(l, m), 2)*(pow(Clm(l, m), 2) + pow(Clm(l - 1, m), 2)) - pow(Clm(l, m), 2)*pow(Clm(l + 1, m), 2);
	}else if(n == 1 && k == 0){
		return l*Clm(l + 1, m);
	}else if(n == 1 && k == 2){
		return Clm(l + 1, m)*(l*(pow(Clm(l + 1, m), 2) + pow(Clm(l + 2, m), 2)) - (l + 1.)*pow(Clm(l, m), 2));
	}else if(n == 2 && k == 1){
		return l*Clm(l + 1, m)*Clm(l + 2, m);
	}else if(n == 2 && k == 3){
		return Clm(l + 1, m)*Clm(l + 2, m)*(l*(pow(Clm(l + 1, m), 2) + pow(Clm(l + 2, m), 2) + pow(Clm(l + 3, m), 2)) - (l + 1.)*pow(Clm(l, m), 2));
	}else if(n == 3 && k == 2){
		return l*Clm(l + 1, m)*Clm(l + 2, m)*Clm(l + 3, m);
	}else if(n == 4 && k == 3){
		return l*Clm(l + 1, m)*Clm(l + 2, m)*Clm(l + 3, m)*Clm(l + 4, m);
	}else{
		return 0.;
	}
}

double xi_coupling_higher_order(int k, double thp){
	double cth = cos(thp);
	double sth = sin(thp);
	switch (k) {
		case 0:
			return 0.5*(2. - 7.*pow(cth, 2) + 8.*pow(cth, 4)*pow(sth, 2))*pow(sth, -7);
		case 1:
			return 1.5*pow(cth, 3)*(1. + 4.*pow(cth, 2))*pow(sth, -7);
		case 2:
			return 0.5*(1. - 8.*pow(cth, 2) - 8.*pow(cth, 4))*pow(sth, -7);
		case 3:
			return 0.5*cth*(3. + 2.*pow(cth, 2))*pow(sth, -7);
		default:
			return 0.;
	}
	return 0.;
}

void beta_coupling_higher_order(Vector &couplingVector, TeukolskyMode &teuk, int n, int l, int m){
	double alpha[4];
	for(int k = 0; k < 4; k++){
		alpha[k] = teuk.getCouplingCoefficient(l)*alpha_coupling_higher_order(n, k, l, m);
	}
	for(size_t jth = 0; jth < couplingVector.size(); jth++){
		for(int k = 0; k < 4; k++){
			couplingVector[jth] += alpha[k]*xi_coupling_higher_order(k, teuk.getPolarPoints(jth));
		}
	}
}

double beta_coupling_higher_order(TeukolskyMode &teuk, int n, int l, int m, int jth){
	double alpha[4];
	for(int k = 0; k < 4; k++){
		alpha[k] = alpha_coupling_higher_order(n, k, l, m);
	}
	double beta = 0.;
	for(int k = 0; k < 4; k++){
		beta += alpha[k]*xi_coupling_higher_order(k, teuk.getPolarPoints(jth));
	}
	return teuk.getCouplingCoefficient(l)*beta;
}

double beta_coupling_higher_order(int n, int l, int m, double theta){
	double alpha[4];
	for(int k = 0; k < 4; k++){
		alpha[k] = alpha_coupling_higher_order(n, k, l, m);
	}
	double beta = 0.;
	for(int k = 0; k < 4; k++){
		beta += alpha[k]*xi_coupling_higher_order(k, theta);
	}
	return beta;
}

void polar_coupling_coefficient_4(Vector &couplingVector, TeukolskyMode &teuk, int l, int m){
	for(int n = -4; n <= 4; n++){
		if(l - n >= abs(m)){
			beta_coupling_higher_order(couplingVector, teuk, n, l - n, m);
		}
	}
}

double polar_coupling_coefficient_4(TeukolskyMode &teuk, int l, int m, int jth){
	double coupling = 0.;
	for(int n = -4; n <= 4; n++){
		if(l - n >= abs(m)){
			coupling += beta_coupling_higher_order(teuk, n, l - n, m, jth);
		}
	}
	return coupling;
}

void polar_coupling_coefficient(int order, Vector &couplingVector, TeukolskyMode &teuk, int l, int m){
	switch(order){
		case 2:
			return polar_coupling_coefficient_2(couplingVector, teuk, l, m);
		case 3:
			return polar_coupling_coefficient_3(couplingVector, teuk, l, m);
		case 4:
			return polar_coupling_coefficient_4(couplingVector, teuk, l, m);
	}
}

double polar_coupling_coefficient(int order, TeukolskyMode &teuk, int l, int m, int jth){
	switch(order){
		case 2:
			return polar_coupling_coefficient_2(teuk, l, m, jth);
		case 3:
			return polar_coupling_coefficient_3(teuk, l, m, jth);
		case 4:
			return polar_coupling_coefficient_4(teuk, l, m, jth);
	}
	return polar_coupling_coefficient_2(teuk, l, m, jth);
}

int radial_n_mode_max_ssf(int l, double){
	return l;
}

RealTensor scalar_self_force_reference(List components, int lmax, int m, GeodesicSource &geo, int sampleSize){
	Vector ref0 = scalar_self_force_reference(components[0], geo, sampleSize);
  RealTensor ssfRef = scalar_self_force_components_data_init(components.size(), lmax, m, ref0.size());
  for(size_t mu = 0; mu < ssfRef.size(); mu++){
    for(size_t l = 0; l < ssfRef[mu].size(); l++){
      ssfRef[mu][l] = scalar_self_force_reference(components[mu], geo, sampleSize);
    }
  }

  return ssfRef;
};

Vector scalar_self_force_reference(int mu, GeodesicSource &geo, int sampleSize){
  if(mu == 0){
    Vector ssfRef = regularization_parameter_Bt_from_source(geo, sampleSize);
    for(size_t i = 0; i < ssfRef.size(); i++){
      if(ssfRef[i] == 0.){
        ssfRef[i] = 1.;
      }
    }
    return ssfRef;
  }else if(mu == 1){
    Vector ssfRef = regularization_parameter_Br_from_source(geo, sampleSize);
    for(size_t i = 0; i < ssfRef.size(); i++){
      if(ssfRef[i] == 0.){
        ssfRef[i] = 1.;
      }
    }
    return ssfRef;
  }else if(mu == 2){
    Vector ssfRef = regularization_parameter_Btheta_from_source(geo, sampleSize);
    for(size_t i = 0; i < ssfRef.size(); i++){
      if(ssfRef[i] == 0.){
        ssfRef[i] = 1.;
      }
    }
    return ssfRef;
  }else if(mu == 3){
    Vector ssfRef = regularization_parameter_Bphi_from_source(geo, sampleSize);
    for(size_t i = 0; i < ssfRef.size(); i++){
      if(ssfRef[i] == 0.){
        ssfRef[i] = 1.;
      }
    }
    return ssfRef;
  }else if(mu == 4){
    Vector ssfRef = regularization_parameter_Bphi_from_source(geo, sampleSize);
    for(size_t i = 0; i < ssfRef.size(); i++){
      if(ssfRef[i] == 0.){
        ssfRef[i] = 1.;
      }
    }
    return ssfRef;
  }else{
    return Vector(sampleSize*sampleSize, 0.);
  }
}

RealTensor scalar_self_force_components_data_init(int componentNum, int lmax, int m, int sampleSize){
  return RealTensor(componentNum, RealMatrix(lmax - abs(m) + 1, Vector(sampleSize, 0.)));
}

RealTensor scalar_self_force_components_convergence_init(const RealTensor &convergenceData, int convergenceCriteria){
  RealTensor convergenceDataCopy(convergenceData.size(), RealMatrix(convergenceData[0].size(), Vector(convergenceData[0][0].size(), 0.)));
  for(size_t mu = 0; mu < convergenceData.size(); mu++){
    for(size_t l = 0; l < convergenceData[0].size(); l++){
      for(size_t i = 0; i < convergenceData[0][0].size(); i++){
        if(convergenceData[mu][l][i] >= convergenceCriteria){
          convergenceDataCopy[mu][l][i] = CONVERGENCE_CRITERIA_LARGE;
        }
      }
    }
  }

  return convergenceDataCopy;
}

SelfForceData scalar_self_force_components_m(int mu, int lmax, int m, GeodesicSource &geo, int sampleSize){
  List components(1);
  components[0] = mu;
  return scalar_self_force_components_m(components, lmax, m, geo, sampleSize);
}

SelfForceData scalar_self_force_components_m(List components, int lmax, int m, GeodesicSource &geo, int sampleSize){
	int fullSampleSize;
	if(geo.getEccentricity() > 0. && abs(geo.getInclination()) < 1.){
		fullSampleSize = sampleSize*sampleSize;
	}else if(geo.getEccentricity() == 0. && abs(geo.getInclination()) == 1.){
		fullSampleSize = 1;
	}else{
		fullSampleSize = sampleSize;
	}
  RealTensor ssfIn = scalar_self_force_components_data_init(components.size(), lmax, m, fullSampleSize);
  RealTensor ssfUp = ssfIn;
  RealTensor ssfRef = scalar_self_force_reference(components, lmax, m, geo, sampleSize);
  scalar_self_force_components_m(components, ssfIn, ssfUp, ssfRef, m, geo);
  SelfForceData ssfData = {
    .in = ssfIn,
    .up = ssfUp,
  };
  return ssfData;
}

int scalar_self_force_components_m(List components, RealTensor &ssfTotIn, RealTensor &ssfTotUp, RealTensor &ssfRef, int m, GeodesicSource geo){

  int componentNum = components.size();
  int lmax = ssfTotIn[0].size() + abs(m) - 1;
  int sampleSize = ssfTotIn[0][0].size();
  int convergenceCriteria = 3;
  convergenceCriteria *= 2; // double it due to the fact that l + m + k = odd modes vanish

  RealTensor ssfInit = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfConvergenceIn = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfConvergenceUp = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfModeIn = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfModeUp = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfModePreviousIn = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfModePreviousUp = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);

  if(abs(geo.getInclination()) == 1.){
		// std::cout << "(SSF) Only calculating k = 0 mode. \n";
    return scalar_self_force_components_mk(components, ssfTotIn, ssfTotUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, ssfRef, m, 0, geo);
  }

  int kInit = 4;
  int kInitMin = -m - kInit;
  int kInitMax = -m + kInit;
  // if(kInitMin > 0){ kInitMin = 0; }
  // if(kInitMax < 0){ kInitMax = 0; }
  int k = kInitMin;
  int kPlusMinus = lmax + 50;

  int check = scalar_self_force_components_mk(components, ssfModeIn, ssfModeUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, ssfRef, m, k, geo);
  int convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
  convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
  k++;
  RealTensor ssfModeInInit = ssfModePreviousIn;
  RealTensor ssfModeUpInit = ssfModePreviousUp;

  // int ltestprint = 20;
  // int kTestIter = 18;
  // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k-1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
  //
  // std::cout << "F_up(l = "<<ltestprint<<") = " << ssfTotUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k-1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";

  while(k <= kInitMax){
    ssfModeIn = ssfInit;
    ssfModeUp = ssfInit;
    check = scalar_self_force_components_mk(components, ssfModeIn, ssfModeUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, ssfRef, m, k, geo);
    convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, CONVERGENCE_CRITERIA_LARGE);
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, CONVERGENCE_CRITERIA_LARGE);
    k++;
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
    //
    // std::cout << "F_up(l = "<<ltestprint<<") = " << ssfTotUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
		if(convergenceFlag != 0){
			std::cout << "(SSF) Error: Large convergence criteria of "<< CONVERGENCE_CRITERIA_LARGE <<" met after summing over only "<< k - kInitMin << " k-modes. \n";
		}
  }

  ssfConvergenceIn = ssfInit;
  ssfConvergenceUp = ssfInit;

  while(convergenceFlag == 0 && k <= kInitMax + kPlusMinus){
    ssfModeIn = ssfInit;
    ssfModeUp = ssfInit;
    check = scalar_self_force_components_mk(components, ssfModeIn, ssfModeUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, ssfRef, m, k, geo);
    convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
    k++;
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
    //
    // std::cout << "F_up(l = "<<ltestprint<<") = " << ssfTotUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k-1<<"\n";
  }

  if(k > kInitMax + kPlusMinus){
    std::cout << "(SSF) ERROR: SSF (m,k) = ("<<m<<")-mode reached kmax = "<<kInitMax + kPlusMinus<<" before converging. \n";
  }

  k = kInitMin - 1;
  ssfModeIn = ssfInit;
  ssfModeUp = ssfInit;
  ssfConvergenceIn = ssfInit;
  ssfConvergenceUp = ssfInit;
  ssfModePreviousIn = ssfModeInInit;
  ssfModePreviousUp = ssfModeUpInit;

  check = scalar_self_force_components_mk(components, ssfModeIn, ssfModeUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, ssfRef, m, k, geo);
  convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
  convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
  k--;
  // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k+1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";
  //
  // std::cout << "F_up(l = "<<ltestprint<<") = " << ssfTotUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k+1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";

  while(convergenceFlag == 0 && k >= kInitMin - kPlusMinus){
    ssfModeIn = ssfInit;
    ssfModeUp = ssfInit;
    check = scalar_self_force_components_mk(components, ssfModeIn, ssfModeUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, ssfRef, m, k, geo);
    convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
    k--;
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k+1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";
    //
    // std::cout << "F_up(l = "<<ltestprint<<") = " << ssfTotUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " with sum from k = "<<kInitMin<<" to "<<k+1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceUp[0][ltestprint - abs(m)][kTestIter*sqrt(sampleSize)] << " for k = "<<k+1<<"\n";
  }

  if(k < kInitMin - kPlusMinus){
    std::cout << "(SSF) ERROR: SSF (m,k) = ("<<m<<")-mode reached kmin = "<<kInitMin - kPlusMinus<<" before converging. \n";
  }

  return convergenceFlag;
}

SelfForceData scalar_self_force_components_mk(List components, int lmax, int m, int k, GeodesicSource &geo, int sampleSize){
	int fullSampleSize;
	if(geo.getEccentricity() > 0. && abs(geo.getInclination()) < 1.){
		fullSampleSize = sampleSize*sampleSize;
	}else if(geo.getEccentricity() == 0. && abs(geo.getInclination()) == 1.){
		fullSampleSize = 1;
	}else{
		fullSampleSize = sampleSize;
	}
  RealTensor ssfIn = scalar_self_force_components_data_init(components.size(), lmax, m, fullSampleSize);
  RealTensor ssfUp = scalar_self_force_components_data_init(components.size(), lmax, m, fullSampleSize);
  RealTensor ssfConvergeIn = ssfIn;
  RealTensor ssfConvergeUp = ssfUp;
  RealTensor ssfRef = scalar_self_force_reference(components, lmax, m, geo, sampleSize);
  scalar_self_force_components_mk(components, ssfIn, ssfUp, ssfConvergeIn, ssfConvergeUp, 5, ssfRef, m, k, geo);
  SelfForceData ssfData = {
    .in = ssfIn,
    .up = ssfUp,
  };
  return ssfData;
}

int scalar_self_force_components_mk(List components, RealTensor &ssfTotIn, RealTensor &ssfTotUp, const RealTensor &ssfConvergenceInM, const RealTensor &ssfConvergenceUpM, int convergenceCriteriaM, RealTensor &ssfRef, int m, int k, GeodesicSource &geo){
	std::cout << "(SSF) Computing (m, k) = ("<<m<<", "<<k<<")\n";
  if(abs(geo.getEccentricity()) == 0.){
		// std::cout << "(SSF) Only calculating n = 0 mode. \n";
    scalar_self_force_components_mkn(components, ssfTotIn, ssfTotUp, m, k, 0, geo);
    return 0;
  }
  int componentNum = components.size();
  int lmax = ssfTotIn[0].size() + abs(m) - 1;
  int sampleSize = ssfTotIn[0][0].size();
  int convergenceCriteria = 5;

  int nMin = minimum_radial_harmonic(m, k, geo);
  // int nInitMin = radial_n_mode_max_ssf(abs(m), geo.getEccentricity());
  // int nInitMax = radial_n_mode_max_ssf(lmax, geo.getEccentricity());
  int nInitMin = -8;
  int nInitMax = 8;
  int nWidth = nInitMax - nInitMin + 1;
  if(m < 0){ nInitMin = -nInitMax; nInitMax = nInitMin + nWidth; }
  if(nInitMin <= nMin){ nInitMin = nMin + 1; nInitMax = nInitMin + nWidth; }
  int n = nInitMin;
	int nbuffer = 120/pow(1. - pow(geo.getEccentricity(), 2), 1.5);
  int nMaxPlus = nbuffer + abs(k);

  RealTensor ssfInit = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfConvergenceInitIn = scalar_self_force_components_convergence_init(ssfConvergenceInM, convergenceCriteriaM);
  RealTensor ssfConvergenceInitUp = scalar_self_force_components_convergence_init(ssfConvergenceUpM, convergenceCriteriaM);
  // ssfConvergenceInitIn = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  // ssfConvergenceInitUp = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfConvergenceIn = ssfConvergenceInitIn;
  RealTensor ssfConvergenceUp = ssfConvergenceInitUp;
  RealTensor ssfModeIn = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfModeUp = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfModePreviousIn = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);
  RealTensor ssfModePreviousUp = scalar_self_force_components_data_init(componentNum, lmax, m, sampleSize);

  int check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, m, k, n, geo);
  if(check != 0){
    std::cout << "(SSF) ERROR: NaN encountered when calculating ("<<m<<","<<k<<","<<n<<") mode \n";
    // return 0;
  }
  int convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
  convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
  n++;
  RealTensor ssfModeInInit = ssfModePreviousIn;
  RealTensor ssfModeUpInit = ssfModePreviousUp;
  // int ltestprint = 20;
  // int iterTest = 18;

  // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n-1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
  //
  // std::cout << "F_up(l = "<<ltestprint<<") = " << ssfTotUp[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n-1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeUp[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceUp[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";

  while(n <= nInitMax){
    ssfModeIn = ssfInit;
    ssfModeUp = ssfInit;
    check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, m, k, n, geo);
    if(check != 0){
      std::cout << "(SSF) ERROR: NaN encountered when calculating ("<<m<<","<<k<<","<<n<<") mode \n";
      // return 0;
    }
    convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, CONVERGENCE_CRITERIA_LARGE);
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, CONVERGENCE_CRITERIA_LARGE);
    n++;
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
    //
    // std::cout << "F_up(l = "<<ltestprint<<") = " << ssfTotUp[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeUp[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceUp[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
		if(convergenceFlag != 0){
			std::cout  << "(SSF) Error: Large convergence criteria of "<< CONVERGENCE_CRITERIA_LARGE <<" met after summing over only "<< n - nInitMin << " n-modes. \n";
		}
  }

  ssfConvergenceIn = ssfConvergenceInitIn;
  ssfConvergenceUp = ssfConvergenceInitUp;
  while(convergenceFlag == 0 && n <= nInitMax + nMaxPlus){
    ssfModeIn = ssfInit;
    ssfModeUp = ssfInit;
    check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, m, k, n, geo);
		// check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, m, k, n, geo);
    if(check != 0){
      std::cout << "(SSF) ERROR: NaN encountered when calculating ("<<m<<","<<k<<","<<n<<") mode \n";
      // return 0;
    }
    convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
    n++;
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
    //
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n-1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n-1<<"\n";
  }
  if(n >= nInitMax + nMaxPlus){
    std::cout << "(SSF) ERROR: SSF (m,k) = ("<<m<<","<<k<<")-mode reached nmax = "<<nInitMax + nMaxPlus<<" before converging. \n";
    for(size_t mu = 0; mu < ssfConvergenceIn.size(); mu++){
      for(size_t i = 0; i < ssfConvergenceIn[mu].size(); i++){
        int lFlag = 0;
        for(size_t j = 0; j < ssfConvergenceIn[mu][i].size(); j++){
          if(ssfConvergenceIn[mu][i][j] < convergenceCriteria){ lFlag += 1; }
          if(ssfConvergenceUp[mu][i][j] < convergenceCriteria){ lFlag += 1; }
        }
        if(lFlag > 0){
          std::cout << "(SSF) ERROR: ("<<i + abs(m)<<","<<m<<","<<k<<")-mode of "<<components[mu]<<"-component did not converge.\n";
        }
      }
    }
  }else{
    std::cout << "(SSF) SSF (m,k) = ("<<m<<", "<<k<<")-mode converged at max n = "<<n-1<<" \n";
  }

  n = nInitMin - 1;
  ssfModeIn = ssfInit;
  ssfModeUp = ssfInit;
  ssfConvergenceIn = ssfConvergenceInitIn;
  ssfConvergenceUp = ssfConvergenceInitUp;
  ssfModePreviousIn = ssfModeInInit;
  ssfModePreviousUp = ssfModeUpInit;
  check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, m, k, n, geo);
  if(check != 0){
    std::cout << "(SSF) ERROR: NaN encountered when calculating ("<<m<<","<<k<<","<<n<<") mode \n";
    // return 0;
  }
  convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
  convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
  n--;
  // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n+1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
  //
  // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n+1<<"\n";
  // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
  // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";

  convergenceFlag = 0;
  if(n >= nMin){
    ssfModeIn = ssfInit;
    ssfModeUp = ssfInit;
    check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, m, k, n, geo);
    if(check != 0){
      std::cout << "(SSF) ERROR: NaN encountered when calculating ("<<m<<","<<k<<","<<n<<") mode \n";
      // return 0;
    }
    convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
    n--;
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n+1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
    //
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n+1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
  }

  while(convergenceFlag == 0 && n >= nMin){
    ssfModeIn = ssfInit;
    ssfModeUp = ssfInit;
    check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, m, k, n, geo);
		// check = scalar_self_force_components_mkn(components, ssfModeIn, ssfModeUp, m, k, n, geo);
    if(check != 0){
      std::cout << "(SSF) ERROR: NaN encountered when calculating ("<<m<<","<<k<<","<<n<<") mode \n";
      // return 0;
    }
    convergenceFlag = scalar_self_force_components_convergence_sum(ssfTotIn, ssfModeIn, ssfModePreviousIn, ssfRef, ssfConvergenceIn, convergenceCriteria);
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTotUp, ssfModeUp, ssfModePreviousUp, ssfRef, ssfConvergenceUp, convergenceCriteria);
    n--;
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n+1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
    //
    // std::cout << "F_in(l = "<<ltestprint<<") = " << ssfTotIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " with sum from n = "<<nInitMin<<" to "<<n+1<<"\n";
    // std::cout << "F_mode(l = "<<ltestprint<<") = " << ssfModeIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
    // std::cout << "convergenceCount = " << ssfConvergenceIn[0][ltestprint - abs(m)][iterTest*sqrt(sampleSize)] << " for n = "<<n+1<<"\n";
  }

  if(n < nMin){
    std::cout << "(SSF) SSF (m,k) = ("<<m<<", "<<k<<")-mode truncated at minimum n = "<<nMin<<" \n";
  }else{
    std::cout << "(SSF) SSF (m,k) = ("<<m<<", "<<k<<")-mode converged at minimum n = "<<n+1<<" \n";
  }

  return convergenceFlag;
}

int scalar_self_force_components_convergence_sum(RealTensor &ssfTot, RealTensor &ssfMode, RealTensor &ssfPrevious, RealTensor &ssfRef, RealTensor &ssfConvergence, int convergenceCriteria){
  int convergenceFlag = 1;

  for(size_t mu = 0; mu < ssfTot.size(); mu++){
    convergenceFlag *= scalar_self_force_components_convergence_sum(ssfTot[mu], ssfMode[mu], ssfPrevious[mu], ssfRef[mu], ssfConvergence[mu], convergenceCriteria);
  }

  return convergenceFlag;
}

int scalar_self_force_components_convergence_sum(RealMatrix &ssfTot, RealMatrix &ssfMode, RealMatrix &ssfPrevious, RealMatrix &ssfRef, RealMatrix &ssfConvergence, int convergenceCriteria){
  int convergenceFlag = 1;

  for(size_t l = 0; l < ssfMode.size(); l++){
    for(size_t i = 0; i < ssfMode[l].size(); i++){
      if(ssfConvergence[l][i] < double(convergenceCriteria)){
        ssfTot[l][i] += ssfMode[l][i];

        if(abs(ssfMode[l][i]/ssfRef[l][i]) < SSF_EPSILON && abs(ssfMode[l][i]) < abs(ssfPrevious[l][i])){
          ssfConvergence[l][i] += 1.;
        }else if(abs(ssfMode[l][i]) < DBL_EPSILON){
          ssfConvergence[l][i] += 1.;
        }else{
          ssfConvergence[l][i] = 0.;
        }

        if(abs(ssfMode[l][i]) > 0.){
          ssfPrevious[l][i] = ssfMode[l][i];
        }

        if(ssfConvergence[l][i] < double(convergenceCriteria)){
          convergenceFlag = 0;
        }
      }
    }
  }

  return convergenceFlag;
}

Complex self_force_amplitude(int mu, TeukolskyMode &teuk, int l, int m, int jr, int jth, BoundaryCondition bc){
  if(mu == 0){
    if( abs(teuk.getFrequency()) > ZERO_FREQ_LIMIT ){
      return -2.*I*teuk.getFrequency()*teuk.getCouplingCoefficient(l)*(teuk.getRadialSolution(bc, jr))*Ylm(l, m, teuk.getPolarPoints(jth));
			// we double count the omega_{mkn} and omega_{-m,-k,-n} except for zero frequencies
    }else{
      return 0.;
    }
  }else if(mu == 1){
    if( abs(teuk.getFrequency()) > ZERO_FREQ_LIMIT ){
      return 2.*teuk.getCouplingCoefficient(l)*teuk.getRadialDerivative(bc, jr)*Ylm(l, m, teuk.getPolarPoints(jth));
    }else{
      return teuk.getCouplingCoefficient(l)*teuk.getRadialDerivative(bc, jr)*Ylm(l, m, teuk.getPolarPoints(jth));
    }
  }else if(mu == 2){
		double polarCoupling = polar_coupling_coefficient(4, teuk, l, m, jth);
    if( abs(teuk.getFrequency()) > ZERO_FREQ_LIMIT ){
      return 2.*polarCoupling*teuk.getRadialSolution(bc, jr)*Ylm(l, m, teuk.getPolarPoints(jth));
    }else{
      return polarCoupling*teuk.getRadialSolution(bc, jr)*Ylm(l, m, teuk.getPolarPoints(jth));
    }
  }else if(mu == 3){
    if(abs(m) > 0){
      return 2.*I*double(m)*teuk.getCouplingCoefficient(l)*(teuk.getRadialSolution(bc, jr))*Ylm(l, m, teuk.getPolarPoints(jth));
    }else{
      return 0.;
    }
  }else if(mu == 4){
    if( abs(teuk.getFrequency()) > ZERO_FREQ_LIMIT ){
      return 2.*teuk.getCouplingCoefficient(l)*(teuk.getRadialSolution(bc, jr));
    }else{
      return teuk.getCouplingCoefficient(l)*(teuk.getRadialSolution(bc, jr));
    }
  }else{
    return 0.;
  }
}

RealMatrix scalar_self_force_components_mkn(int mu, int lmax, int m, int k, int n, GeodesicSource &geo, int sampleSize){
	int fullSampleSize;
	if(geo.getEccentricity() > 0. && abs(geo.getInclination()) < 1.){
		fullSampleSize = sampleSize*sampleSize;
	}else if(geo.getEccentricity() == 0. && abs(geo.getInclination()) == 1.){
		fullSampleSize = 1;
	}else{
		fullSampleSize = sampleSize;
	}
  RealTensor ssf = scalar_self_force_components_data_init(1, lmax, m, fullSampleSize);
  RealTensor ssf2 = ssf;
  List components(1);
  components[0] = mu;
  scalar_self_force_components_mkn(components, ssf, ssf2, m, k, n, geo);
  return ssf[0];
}

SelfForceData scalar_self_force_components_mkn(List components, int lmax, int m, int k, int n, GeodesicSource &geo, int sampleSize){
	int fullSampleSize;
	if(geo.getEccentricity() > 0. && abs(geo.getInclination()) < 1.){
		fullSampleSize = sampleSize*sampleSize;
	}else if(geo.getEccentricity() == 0. && abs(geo.getInclination()) == 1.){
		fullSampleSize = 1;
	}else{
		fullSampleSize = sampleSize;
	}
  RealTensor ssfIn = scalar_self_force_components_data_init(components.size(), lmax, m, fullSampleSize);
  RealTensor ssfUp = ssfIn;
  scalar_self_force_components_mkn(components, ssfIn, ssfUp, m, k, n, geo);
  SelfForceData ssfData = {
    .in = ssfIn,
    .up = ssfUp
  };
  return ssfData;
}

int scalar_self_force_components_mkn(List components, RealTensor &ssfIn, RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo){
	if(geo.getEccentricity() == 0. && abs(geo.getInclination()) == 1.){
		return scalar_self_force_components_mkn_circular(components, ssfIn, ssfUp, m, k, n, geo);
	}else if(abs(geo.getInclination()) == 1.){
		return scalar_self_force_components_mkn_equatorial(components, ssfIn, ssfUp, m, k, n, geo);
	}else if( geo.getEccentricity() == 0. ){
		return scalar_self_force_components_mkn_spherical(components, ssfIn, ssfUp, m, k, n, geo);
	}

	return scalar_self_force_components_mkn_generic(components, ssfIn, ssfUp, m, k, n, geo);
}

int scalar_self_force_components_mkn_equatorial(List components, RealTensor &ssfIn, RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo){
	int returnFlag = 0;
	if(abs(k) > 0){
		return returnFlag;
	}

  int componentNum = components.size();
  int lmax = ssfIn[0].size() - 1 + abs(m);
  double freq = geo.getTimeFrequency(m, k, n);
  SpinWeightedHarmonic swsh(0, lmax, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmax = lmax;
  while(abs(swsh.getCouplingCoefficient()[jmax - abs(m)]) > DBL_EPSILON && jmax < swsh.getMaxCouplingModeNumber() - 1){
    jmax += 2;
  }

  Vector deltaTR = geo.getTimeAccumulation(1);
  Vector deltaPhiR = geo.getAzimuthalAccumulation(1);
  int sampleSize = ssfIn[0][0].size();
  int sampleSizeHalf = sampleSize/2;
  int radialSampleSize = deltaTR.size() - 1;
  int radialSampleRate = radialSampleSize/sampleSizeHalf;

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  double phaseR = 0.;
  for(int j = abs(m); j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      teuk.generateSolutions(geo);
			// if(j < 40){
			//   std::cout << "Teukolsky mode ("<<j<<","<<m<<","<<k<<","<<n<<") computed \n";
			//   std::cout << "ClmUp = " << teuk.getTeukolskyAmplitude(Up) << "\n";
			//   std::cout << "ClmIn = " << teuk.getTeukolskyAmplitude(In)*sqrt(2.*(1. + sqrt(1. - pow(teuk.getBlackHoleSpin(),2)))) << "\n";
			// }
      int couplingMax = teuk.getCouplingCoefficient().size() + abs(m) - 1;
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = abs(m); l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
					int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            calcFlag = 0;
          }
          if(calcFlag == 1){
            for(int jr = 0; jr <= sampleSizeHalf; jr++){
              ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, 0, Up);
              ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, 0, In);

              phaseR = freq*deltaTR[jr*radialSampleRate] - double(m)*deltaPhiR[jr*radialSampleRate] + 2.*M_PI*double(n*jr)/double(sampleSize);

              if( isnan(abs(ssfModeTermIn)) ){
                if(returnFlag < 2){
                  std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                  returnFlag = 2;
                }
              }else{
                ssfIn[mu][l - abs(m)][jr] += std::real(ssfModeTermIn*exp(-I*phaseR));
                if(jr > 0 && jr < sampleSizeHalf){
                  ssfIn[mu][l - abs(m)][(sampleSize - jr)] += std::real(ssfModeTermIn*exp(I*phaseR));
                }
              }

              if( isnan(abs(ssfModeTermUp)) ){
                if(returnFlag < 2){
                  std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                  returnFlag = 2;
                }
              }else{
                ssfUp[mu][l - abs(m)][jr] += std::real(ssfModeTermUp*exp(-I*phaseR));
                if(jr > 0 && jr < sampleSizeHalf){
                  ssfUp[mu][l - abs(m)][(sampleSize - jr)] += std::real(ssfModeTermUp*exp(I*phaseR));
                }
              }
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

int scalar_self_force_components_mkn_spherical(List components, RealTensor &ssfIn, RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo){
	int returnFlag = 0;
	if(abs(n) > 0){
		return returnFlag;
	}

  int componentNum = components.size();
  int lmax = ssfIn[0].size() - 1 + abs(m);
  double freq = geo.getTimeFrequency(m, k, n);
  SpinWeightedHarmonic swsh(0, lmax, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmax = lmax;
  while(abs(swsh.getCouplingCoefficient()[jmax - abs(m)]) > DBL_EPSILON && jmax < swsh.getMaxCouplingModeNumber() - 1){
    jmax += 2;
  }

  Vector deltaTTh = geo.getTimeAccumulation(2);
  Vector deltaPhiTh = geo.getAzimuthalAccumulation(2);
  int sampleSize = ssfIn[0][0].size();
  int sampleSizeHalf = sampleSize/2;
  int polarSampleSize = deltaTTh.size() - 1;
  int polarSampleRate = polarSampleSize/sampleSizeHalf;

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  double phaseTh = 0.;
  for(int j = abs(m); j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      teuk.generateSolutions(geo);
      int couplingMax = teuk.getCouplingCoefficient().size() + abs(m) - 1;
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = abs(m); l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
					int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            // if( (l + j) % 2 == 0 ){
            //   calcFlag = 0;
            // }else if(abs(teuk.getCouplingCoefficient(l - 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l - 1)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 1)) < COUPLING_EPSILON){
            //   calcFlag = 0;
            // }
          }
          if(calcFlag == 1){
            for(int jth = 0; jth <= sampleSizeHalf; jth++){
              ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, 0, jth*polarSampleRate, Up);
              ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, 0, jth*polarSampleRate, In);

            	phaseTh = freq*deltaTTh[jth*polarSampleRate] - double(m)*deltaPhiTh[jth*polarSampleRate] + 2.*M_PI*double(k*jth)/double(sampleSize);

            	if( isnan(abs(ssfModeTermIn)) ){
              	if(returnFlag < 2){
                	std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                	returnFlag = 2;
            		}
            	}else{
              	ssfIn[mu][l - abs(m)][jth] += std::real(ssfModeTermIn*exp(-I*phaseTh));
              	if(jth > 0 && jth < sampleSizeHalf){
                	ssfIn[mu][l - abs(m)][sampleSize - jth] += std::real(ssfModeTermIn*exp(I*phaseTh));
                }
              }

              if( isnan(abs(ssfModeTermUp)) ){
                if(returnFlag < 2){
                  std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                  returnFlag = 2;
                }
              }else{
                ssfUp[mu][l - abs(m)][jth] += std::real(ssfModeTermUp*exp(-I*phaseTh));
                if(jth > 0 && jth < sampleSizeHalf){
                  ssfUp[mu][l - abs(m)][sampleSize - jth] += std::real(ssfModeTermUp*exp(I*phaseTh));
                }
              }
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

int scalar_self_force_components_mkn_circular(List components, RealTensor &ssfIn, RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo){
	int returnFlag = 0;
	if(abs(k) > 0 || abs(n) > 0){
		return returnFlag;
	}

  int componentNum = components.size();
  int lmax = ssfIn[0].size() - 1 + abs(m);
  double freq = geo.getTimeFrequency(m, k, n);
  SpinWeightedHarmonic swsh(0, lmax, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmax = lmax;
  while(abs(swsh.getCouplingCoefficient()[jmax - abs(m)]) > DBL_EPSILON && jmax < swsh.getMaxCouplingModeNumber() - 1){
    jmax += 2;
  }

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  for(int j = abs(m); j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      teuk.generateSolutions(geo);
      int couplingMax = teuk.getCouplingCoefficient().size() + abs(m) - 1;
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = abs(m); l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
					int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            calcFlag = 0;
          }
          if(calcFlag == 1){
						ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, 0, 0, Up);
						ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, 0, 0, In);

            if( isnan(abs(ssfModeTermIn)) ){
              if(returnFlag < 2){
                std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                returnFlag = 2;
              }
            }else{
              ssfIn[mu][l - abs(m)][0] += std::real(ssfModeTermIn);
            }

            if( isnan(abs(ssfModeTermUp)) ){
              if(returnFlag < 2){
                std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                returnFlag = 2;
              }
            }else{
              ssfUp[mu][l - abs(m)][0] += std::real(ssfModeTermUp);
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

int scalar_self_force_components_mkn_generic(List components, RealTensor &ssfIn, RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo){
	int returnFlag = 0;

	if(geo.getEccentricity() == 0.){
		if(abs(n) > 0){
			return returnFlag;
		}else if(abs(geo.getInclination()) == 1.){
			if(abs(k) > 0){
				return returnFlag;
			}else{
				return scalar_self_force_components_mkn_circular(components, ssfIn, ssfUp, m, k, n, geo);
			}
		}else{
			return scalar_self_force_components_mkn_spherical(components, ssfIn, ssfUp, m, k, n, geo);
		}
  }else if(abs(geo.getInclination()) == 1.){
		if(abs(k) > 0){
			return returnFlag;
		}else{
			return scalar_self_force_components_mkn_equatorial(components, ssfIn, ssfUp, m, k, n, geo);
		}
  }

  int componentNum = components.size();
  int lmax = ssfIn[0].size() - 1 + abs(m);
  double freq = geo.getTimeFrequency(m, k, n);
	int lmaxTest = lmax + 4; // we include the plus 4 to account for the extra coupling for the theta component
  SpinWeightedHarmonic swsh(0, lmaxTest, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmax = lmaxTest;
  while(abs(swsh.getCouplingCoefficient()[jmax - abs(m)]) > DBL_EPSILON && jmax < swsh.getMaxCouplingModeNumber() - 1){
    jmax += 2;
  }

  Vector deltaTR = geo.getTimeAccumulation(1);
  Vector deltaTTh = geo.getTimeAccumulation(2);
  Vector deltaPhiR = geo.getAzimuthalAccumulation(1);
  Vector deltaPhiTh = geo.getAzimuthalAccumulation(2);
  int sampleSize = sqrt(ssfIn[0][0].size());
  int sampleSizeHalf = sampleSize/2;
  int radialSampleSize = deltaTR.size() - 1;
  int radialSampleRate = radialSampleSize/sampleSizeHalf;
  int polarSampleSize = deltaTTh.size() - 1;
  int polarSampleRate = polarSampleSize/sampleSizeHalf;

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  double phaseR = 0.;
  double phaseTh = 0.;
  for(int j = abs(m); j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      teuk.generateSolutions(geo);
      // if(j < 40){
      //   std::cout << "Teukolsky mode ("<<j<<","<<m<<","<<k<<","<<n<<") computed \n";
      //   std::cout << "ClmUp = " << teuk.getTeukolskyAmplitude(Up) << "\n";
      //   std::cout << "ClmIn = " << teuk.getTeukolskyAmplitude(In)*sqrt(2.*(1. + sqrt(1. - pow(teuk.getBlackHoleSpin(),2)))) << "\n";
      // }
      int couplingMax = teuk.getCouplingCoefficient().size() + abs(m) - 1;
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = abs(m); l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
					int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            // if( (l + j) % 2 == 0 ){
            //   calcFlag = 0;
            // }else if(abs(teuk.getCouplingCoefficient(l - 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l - 1)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 1)) < COUPLING_EPSILON){
            //   calcFlag = 0;
            // }
          }
          if(calcFlag == 1){
            for(int jr = 0; jr <= sampleSizeHalf; jr++){
              for(int jth = 0; jth <= sampleSizeHalf; jth++){
                ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, jth*polarSampleRate, Up);
                ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, jth*polarSampleRate, In);

                phaseR = freq*deltaTR[jr*radialSampleRate] - double(m)*deltaPhiR[jr*radialSampleRate] + 2.*M_PI*double(n*jr)/double(sampleSize);
                phaseTh = freq*deltaTTh[jth*polarSampleRate] - double(m)*deltaPhiTh[jth*polarSampleRate] + 2.*M_PI*double(k*jth)/double(sampleSize);

                if( isnan(abs(ssfModeTermIn)) ){
                  if(returnFlag < 2){
                    std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                    returnFlag = 2;
                  }
                }else{
                  ssfIn[mu][l - abs(m)][jr*sampleSize + jth] += std::real(ssfModeTermIn*exp(-I*(phaseR + phaseTh)));
                  if(jr > 0 && jth > 0 && jr < sampleSizeHalf && jth < sampleSizeHalf){
                    ssfIn[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermIn*exp(-I*(phaseR - phaseTh)));
                    ssfIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermIn*exp(I*(phaseR - phaseTh)));
                    ssfIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + sampleSize - jth] += std::real(ssfModeTermIn*exp(I*(phaseR + phaseTh)));
                  }else if(jr > 0 && jr < sampleSizeHalf){
                    ssfIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermIn*exp(I*(phaseR - phaseTh)));
                  }else if(jth > 0 && jth < sampleSizeHalf){
                    ssfIn[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermIn*exp(-I*(phaseR - phaseTh)));
                  }
                }

                if( isnan(abs(ssfModeTermUp)) ){
                  if(returnFlag < 2){
                    std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                    returnFlag = 2;
                  }
                }else{
                  ssfUp[mu][l - abs(m)][jr*sampleSize + jth] += std::real(ssfModeTermUp*exp(-I*(phaseR + phaseTh)));
                  if(jr > 0 && jth > 0 && jr < sampleSizeHalf && jth < sampleSizeHalf){
                    ssfUp[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermUp*exp(-I*(phaseR - phaseTh)));
                    ssfUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermUp*exp(I*(phaseR - phaseTh)));
                    ssfUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + sampleSize - jth] += std::real(ssfModeTermUp*exp(I*(phaseR + phaseTh)));
                  }else if(jr > 0 && jr < sampleSizeHalf){
                    ssfUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermUp*exp(I*(phaseR - phaseTh)));
                  }else if(jth > 0 && jth < sampleSizeHalf){
                    ssfUp[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermUp*exp(-I*(phaseR - phaseTh)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

int min_mode_convergence(const RealTensor &ssfConvergence, int convergenceCriteria){
  int minL = -1;
  int lmax = ssfConvergence[0].size();
  int l = 0;
  while(l < lmax){
    int converge = 1;
    for(size_t j = 0; j < ssfConvergence[0][0].size(); j++){
      for(size_t mu = 0; mu < ssfConvergence.size(); mu++){
        if(ssfConvergence[mu][l][j] < convergenceCriteria){
          converge = 0;
        }
      }
    }
    if(converge > 0 && l == minL + 1){
      minL = l;
    }else{
      l = lmax + 1;
    }
    l++;
  }

  return minL + 1;
}

int max_mode_convergence(const RealTensor &ssfConvergence, int convergenceCriteria){
  int maxL = ssfConvergence[0].size();
  int l = maxL - 1;
  while(l > 0){
    int converge = 1;
    for(size_t j = 0; j < ssfConvergence[0][0].size(); j++){
      for(size_t mu = 0; mu < ssfConvergence.size(); mu++){
        if(ssfConvergence[mu][l][j] < convergenceCriteria){
          converge = 0;
        }
      }
    }
    if(converge > 0 && l == maxL - 1){
      maxL = l;
    }else{
      l = -1;
    }
    l--;
  }

  return maxL - 1;
}

int tensor_same_size(const RealTensor &ssf1, const RealTensor &ssf2){
  if(ssf1.size() != ssf2.size()){
    return 0;
  }else if(ssf1[0].size() != ssf2[0].size()){
    return 0;
  }else if(ssf1[0][0].size() != ssf2[0][0].size()){
    return 0;
  }else{
    return 1;
  }
}

int scalar_self_force_components_mkn(List components, RealTensor &ssfIn, RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp, int convergenceCriteria, int m, int k, int n, GeodesicSource &geo){
	if(geo.getEccentricity() == 0. && abs(geo.getInclination()) == 1.){
		return scalar_self_force_components_mkn_circular(components, ssfIn, ssfUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, m, k, n, geo);
	}else if(abs(geo.getInclination()) == 1.){
		return scalar_self_force_components_mkn_equatorial(components, ssfIn, ssfUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, m, k, n, geo);
	}else if( geo.getEccentricity() == 0. ){
		return scalar_self_force_components_mkn_spherical(components, ssfIn, ssfUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, m, k, n, geo);
	}

	return scalar_self_force_components_mkn_generic(components, ssfIn, ssfUp, ssfConvergenceIn, ssfConvergenceUp, convergenceCriteria, m, k, n, geo);
}

int scalar_self_force_components_mkn_equatorial(List components, RealTensor &ssfIn, RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp, int convergenceCriteria, int m, int k, int n, GeodesicSource &geo){
  int returnFlag = 0;
	if(abs(k) > 0){
		return returnFlag;
	}

  double freq = geo.getTimeFrequency(m, k, n);
  int componentNum = components.size();

  int lmaxIn = max_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lmaxUp = max_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmax = lmaxIn < lmaxUp ? lmaxUp : lmaxIn;
  lmax += abs(m);
	// int lmax = ssfIn[0].size() - 1 + abs(m);

  SpinWeightedHarmonic swsh(0, lmax, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmax = lmax;
  while(jmax < swsh.getMaxCouplingModeNumber() - 1 && abs(swsh.getCouplingCoefficient(jmax)) > DBL_EPSILON){
    jmax += 2;
  }

  int lminIn = min_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lminUp = min_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmin = lminIn < lminUp ? lminIn : lminUp;
  lmin += abs(m);
	// int lmin = abs(m);
  swsh = SpinWeightedHarmonic(0, lmin, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmin = lmin;
  while(jmin > swsh.getMinCouplingModeNumber() + 1 && abs(swsh.getCouplingCoefficient(jmin)) > DBL_EPSILON){
    jmin -= 2;
  }

  Vector deltaTR = geo.getTimeAccumulation(1);
  Vector deltaPhiR = geo.getAzimuthalAccumulation(1);
	int sampleSize = ssfIn[0][0].size();
  int sampleSizeHalf = sampleSize/2;
  int radialSampleSize = deltaTR.size() - 1;
  int radialSampleRate = radialSampleSize/sampleSizeHalf;

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  double phaseR = 0.;
  double phaseTh = 0.;
  for(int j = jmin; j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      teuk.generateSolutions(geo);
      int couplingMax = teuk.getMaxCouplingModeNumber();
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = lmin; l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
          int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            calcFlag = 0;
          }
          if(calcFlag == 1){
            for(int jr = 0; jr <= sampleSizeHalf; jr++){
              phaseR = freq*deltaTR[jr*radialSampleRate] - double(m)*deltaPhiR[jr*radialSampleRate] + 2.*M_PI*double(n*jr)/double(sampleSize);

							int calcFlagIn = 0;
							if(abs(teuk.getTeukolskyAmplitude(In)) > 0.){
								if(ssfConvergenceIn[mu][l - abs(m)][jr] < double(convergenceCriteria)){
									calcFlagIn = 1;
								}

								if(jr > 0 && jr < sampleSizeHalf){
									if(ssfConvergenceIn[mu][l - abs(m)][(sampleSize - jr)] < double(convergenceCriteria)){
										calcFlagIn = 1;
									}
								}
							}

              if(calcFlagIn == 1){
                ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, 0, In);

                if( isnan(abs(ssfModeTermIn)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                  if(returnFlag < 2){
                    std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                    returnFlag = 2;
                  }
                }else{
                  ssfIn[mu][l - abs(m)][jr] += std::real(ssfModeTermIn*exp(-I*(phaseR)));
                  if(jr > 0 && jr < sampleSizeHalf){
                    ssfIn[mu][l - abs(m)][(sampleSize - jr)] += std::real(ssfModeTermIn*exp(I*(phaseR)));
                  }
                }
              }

							int calcFlagUp = 0;
							if(abs(teuk.getTeukolskyAmplitude(Up)) > 0.){
								if(ssfConvergenceUp[mu][l - abs(m)][jr] < double(convergenceCriteria)){
									calcFlagUp = 1;
								}

								if(jr > 0 && jr < sampleSizeHalf){
									if(ssfConvergenceUp[mu][l - abs(m)][(sampleSize - jr)] < double(convergenceCriteria)){
										calcFlagUp = 1;
									}
								}
							}

              if(calcFlagUp == 1){
                ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, 0, Up);

                if( isnan(abs(ssfModeTermUp)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                  if(returnFlag < 2){
                    std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                    returnFlag = 2;
                  }
                }else{
                  ssfUp[mu][l - abs(m)][jr] += std::real(ssfModeTermUp*exp(-I*(phaseR)));
                  if(jr > 0 && jr < sampleSizeHalf){
                    ssfUp[mu][l - abs(m)][(sampleSize - jr)] += std::real(ssfModeTermUp*exp(I*(phaseR)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

int scalar_self_force_components_mkn_spherical(List components, RealTensor &ssfIn, RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp, int convergenceCriteria, int m, int k, int n, GeodesicSource &geo){
  int returnFlag = 0;
	if(abs(n) > 0){
		return returnFlag;
	}

  double freq = geo.getTimeFrequency(m, k, n);
  int componentNum = components.size();

  int lmaxIn = max_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lmaxUp = max_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmax = lmaxIn < lmaxUp ? lmaxUp : lmaxIn;
  lmax += abs(m);
	// int lmax = ssfIn[0].size() - 1 + abs(m);

	int lmaxTest = lmax + 4; // we include the plus 4 to account for the extra coupling for the theta component
	SpinWeightedHarmonic swsh(0, lmaxTest, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
	swsh.generateCouplingCoefficients();
	int jmax = lmaxTest;
	while(jmax < swsh.getMaxCouplingModeNumber() - 1 && abs(swsh.getCouplingCoefficient(jmax)) > DBL_EPSILON){
    jmax += 2;
  }

  int lminIn = min_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lminUp = min_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmin = lminIn < lminUp ? lminIn : lminUp;
  lmin += abs(m);
	// int lmin = abs(m);
	int lminTest = lmin - 4; // we include the minus 4 to account for the extra coupling for the theta component
	int jmin = lminTest;
	if(lminTest < abs(m)){
		lminTest = abs(m);
		jmin = lminTest;
	}else{
		swsh = SpinWeightedHarmonic(0, lminTest, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
	  swsh.generateCouplingCoefficients();
	  while(jmin > swsh.getMinCouplingModeNumber() + 1 && abs(swsh.getCouplingCoefficient(jmin)) > DBL_EPSILON){
	    jmin -= 2;
	  }
	}

  Vector deltaTTh = geo.getTimeAccumulation(2);
  Vector deltaPhiTh = geo.getAzimuthalAccumulation(2);
	int sampleSize = ssfIn[0][0].size();
  int sampleSizeHalf = sampleSize/2;
  int polarSampleSize = deltaTTh.size() - 1;
  int polarSampleRate = polarSampleSize/sampleSizeHalf;

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  double phaseR = 0.;
  double phaseTh = 0.;
  for(int j = jmin; j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      teuk.generateSolutions(geo);
      int couplingMax = teuk.getMaxCouplingModeNumber();
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = lmin; l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
          int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            // if( (l + j) % 2 == 0 ){
            //   calcFlag = 0;
            // }else if(abs(teuk.getCouplingCoefficient(l - 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l - 1)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 1)) < COUPLING_EPSILON){
            //   calcFlag = 0;
            // }
          }
          if(calcFlag == 1){
            for(int jth = 0; jth <= sampleSizeHalf; jth++){
              phaseTh = freq*deltaTTh[jth*polarSampleRate] - double(m)*deltaPhiTh[jth*polarSampleRate] + 2.*M_PI*double(k*jth)/double(sampleSize);

							int calcFlagIn = 0;
							if(abs(teuk.getTeukolskyAmplitude(In)) > 0.){
								if(ssfConvergenceIn[mu][l - abs(m)][jth] < double(convergenceCriteria)){
									calcFlagIn = 1;
								}

								if(jth > 0 && jth < sampleSizeHalf){
									if(ssfConvergenceIn[mu][l - abs(m)][sampleSize - jth] < double(convergenceCriteria)){
										calcFlagIn = 1;
									}
								}
							}

              if(calcFlagIn == 1){
                ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, 0, jth*polarSampleRate, In);

                if( isnan(abs(ssfModeTermIn)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                  if(returnFlag < 2){
                    std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                    returnFlag = 2;
                  }
                }else{
                  ssfIn[mu][l - abs(m)][jth] += std::real(ssfModeTermIn*exp(-I*(phaseTh)));
                  if(jth > 0 && jth < sampleSizeHalf){
                    ssfIn[mu][l - abs(m)][sampleSize - jth] += std::real(ssfModeTermIn*exp(I*(phaseTh)));
                  }
                }
              }

							int calcFlagUp = 0;
							if(abs(teuk.getTeukolskyAmplitude(Up)) > 0.){
								if(ssfConvergenceUp[mu][l - abs(m)][jth] < double(convergenceCriteria)){
									calcFlagUp = 1;
								}

								if(jth > 0 && jth < sampleSizeHalf){
									if(ssfConvergenceUp[mu][l - abs(m)][sampleSize - jth] < double(convergenceCriteria)){
										calcFlagUp = 1;
									}
								}
							}

              if(calcFlagUp == 1){
                ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, 0, jth*polarSampleRate, Up);

                if( isnan(abs(ssfModeTermUp)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                  if(returnFlag < 2){
                    std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                    returnFlag = 2;
                  }
                }else{
                  ssfUp[mu][l - abs(m)][jth] += std::real(ssfModeTermUp*exp(-I*(phaseTh)));
                  if(jth > 0 && jth < sampleSizeHalf){
                    ssfUp[mu][l - abs(m)][sampleSize - jth] += std::real(ssfModeTermUp*exp(I*(phaseTh)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

int scalar_self_force_components_mkn_circular(List components, RealTensor &ssfIn, RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp, int convergenceCriteria, int m, int k, int n, GeodesicSource &geo){
  int returnFlag = 0;
	if(abs(k) > 0 || abs(n) > 0){
		return returnFlag;
	}

  double freq = geo.getTimeFrequency(m, k, n);
  int componentNum = components.size();

  int lmaxIn = max_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lmaxUp = max_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmax = lmaxIn < lmaxUp ? lmaxUp : lmaxIn;
  lmax += abs(m);
	// int lmax = ssfIn[0].size() - 1 + abs(m);

  SpinWeightedHarmonic swsh(0, lmax, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmax = lmax;
  while(jmax < swsh.getMaxCouplingModeNumber() - 1 && abs(swsh.getCouplingCoefficient(jmax)) > DBL_EPSILON){
    jmax += 2;
  }

  int lminIn = min_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lminUp = min_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmin = lminIn < lminUp ? lminIn : lminUp;
  lmin += abs(m);
	// int lmin = abs(m);
  swsh = SpinWeightedHarmonic(0, lmin, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
  swsh.generateCouplingCoefficients();
  int jmin = lmin;
  while(jmin > swsh.getMinCouplingModeNumber() + 1 && abs(swsh.getCouplingCoefficient(jmin)) > DBL_EPSILON){
    jmin -= 2;
  }

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  double phaseR = 0.;
  double phaseTh = 0.;
  for(int j = jmin; j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      teuk.generateSolutions(geo);
      int couplingMax = teuk.getMaxCouplingModeNumber();
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = lmin; l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
          int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            calcFlag = 0;
          }
          if(calcFlag == 1){
						int calcFlagIn = 0;
						if(abs(teuk.getTeukolskyAmplitude(In)) > 0.){
							if(ssfConvergenceIn[mu][l - abs(m)][0] < double(convergenceCriteria)){
								calcFlagIn = 1;
							}
						}

            if(calcFlagIn == 1){
              ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, 0, 0, In);

              if( isnan(abs(ssfModeTermIn)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                if(returnFlag < 2){
                  std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                  returnFlag = 2;
                }
              }else{
                ssfIn[mu][l - abs(m)][0] += std::real(ssfModeTermIn);
              }
            }

						int calcFlagUp = 0;
						if(abs(teuk.getTeukolskyAmplitude(Up)) > 0.){
							if(ssfConvergenceUp[mu][l - abs(m)][0] < double(convergenceCriteria)){
								calcFlagUp = 1;
							}
						}

            if(calcFlagUp == 1){
              ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, 0, 0, Up);

              if( isnan(abs(ssfModeTermUp)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                if(returnFlag < 2){
                  std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                  returnFlag = 2;
                }
              }else{
                ssfUp[mu][l - abs(m)][0] += std::real(ssfModeTermUp);
              }
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

int scalar_self_force_components_mkn_generic(List components, RealTensor &ssfIn, RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp, int convergenceCriteria, int m, int k, int n, GeodesicSource &geo){
  if(geo.getEccentricity() == 0 && abs(n) > 0){
    return 0;
  }else if(abs(geo.getInclination()) == 1 && abs(k) > 0){
    return 0;
  }
  int returnFlag = 0;

  double freq = geo.getTimeFrequency(m, k, n);
  int componentNum = components.size();

  int lmaxIn = max_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lmaxUp = max_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmax = lmaxIn < lmaxUp ? lmaxUp : lmaxIn;
  lmax += abs(m);
	// int lmax = ssfIn[0].size() - 1 + abs(m);

	int lmaxTest = lmax + 4; // we include the plus 4 to account for the extra coupling for the theta component
	SpinWeightedHarmonic swsh(0, lmaxTest, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
	swsh.generateCouplingCoefficients();
	int jmax = lmaxTest;
	while(jmax < swsh.getMaxCouplingModeNumber() - 1 && abs(swsh.getCouplingCoefficient(jmax)) > DBL_EPSILON){
    jmax += 2;
  }

  int lminIn = min_mode_convergence(ssfConvergenceIn, convergenceCriteria);
  int lminUp = min_mode_convergence(ssfConvergenceUp, convergenceCriteria);
  int lmin = lminIn < lminUp ? lminIn : lminUp;
  lmin += abs(m);
	// int lmin = abs(m);
	int lminTest = lmin - 4; // we include the minus 4 to account for the extra coupling for the theta component
	int jmin = lminTest;
	if(lminTest < abs(m)){
		lminTest = abs(m);
		jmin = lminTest;
	}else{
		swsh = SpinWeightedHarmonic(0, lminTest, m, freq*geo.getBlackHoleSpin(), geo.getPolarPosition());
	  swsh.generateCouplingCoefficients();
	  while(jmin > swsh.getMinCouplingModeNumber() + 1 && abs(swsh.getCouplingCoefficient(jmin)) > DBL_EPSILON){
	    jmin -= 2;
	  }
	}

  Vector deltaTR = geo.getTimeAccumulation(1);
  Vector deltaTTh = geo.getTimeAccumulation(2);
  Vector deltaPhiR = geo.getAzimuthalAccumulation(1);
  Vector deltaPhiTh = geo.getAzimuthalAccumulation(2);
	int sampleSizeSquared = ssfIn[0][0].size();
  int sampleSize = sqrt(sampleSizeSquared);
	if(sampleSize*sampleSize != sampleSizeSquared){
		std::cout << "(SSF) Error in the sampling of the self-force structures \n";
	}
  int sampleSizeHalf = sampleSize/2;
  int radialSampleSize = deltaTR.size() - 1;
  int radialSampleRate = radialSampleSize/sampleSizeHalf;
  int polarSampleSize = deltaTTh.size() - 1;
  int polarSampleRate = polarSampleSize/sampleSizeHalf;

  Complex ssfModeTermUp = 0., ssfModeTermIn = 0.;
  double phaseR = 0.;
  double phaseTh = 0.;
  for(int j = jmin; j <= jmax; j++){
    if((abs(m) + j + k) % 2 == 0){
      TeukolskyMode teuk(0, j, m, k, n, geo);
      // std::cout << "Computing Teukolsky mode ("<<j<<","<<m<<","<<k<<","<<n<<")\n";
      teuk.generateSolutions(geo);
      // if(j < 40){
      //   std::cout << "ClmUp = " << teuk.getTeukolskyAmplitude(Up) << "\n";
      //   std::cout << "ClmIn = " << teuk.getTeukolskyAmplitude(In)*sqrt(2.*(1. + sqrt(1. - pow(teuk.getBlackHoleSpin(),2)))) << "\n";
      // }
      int couplingMax = teuk.getMaxCouplingModeNumber();
      int lmaxJMode = lmax < couplingMax ? lmax : couplingMax;
      for(int l = lmin; l <= lmaxJMode; l++){
        for(int mu = 0; mu < componentNum; mu++){
          int calcFlag = 1;
          if(abs(teuk.getTeukolskyAmplitude(Up)) + abs(teuk.getTeukolskyAmplitude(In)) == 0.){
            calcFlag = 0;
          }else if(components[mu] != 2){
            if((l + j) % 2 == 1 || abs(teuk.getCouplingCoefficient(l)) < COUPLING_EPSILON){
              calcFlag = 0;
            }else if(components[mu] == 3 && abs(m) == 0){
              calcFlag = 0;
            }
          }else if(components[mu] == 2){
            // if( (l + j) % 2 == 0 ){
            //   calcFlag = 0;
            // }else if(abs(teuk.getCouplingCoefficient(l - 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l - 1)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 3)) < COUPLING_EPSILON && abs(teuk.getCouplingCoefficient(l + 1)) < COUPLING_EPSILON){
            //   calcFlag = 0;
            // }
          }
          if(calcFlag == 1){
            for(int jr = 0; jr <= sampleSizeHalf; jr++){
              for(int jth = 0; jth <= sampleSizeHalf; jth++){
                phaseR = freq*deltaTR[jr*radialSampleRate] - double(m)*deltaPhiR[jr*radialSampleRate] + 2.*M_PI*double(n*jr)/double(sampleSize);
                phaseTh = freq*deltaTTh[jth*polarSampleRate] - double(m)*deltaPhiTh[jth*polarSampleRate] + 2.*M_PI*double(k*jth)/double(sampleSize);

								int calcFlagIn = 0;
								if(abs(teuk.getTeukolskyAmplitude(In)) > 0.){
									if(ssfConvergenceIn[mu][l - abs(m)][jr*sampleSize + jth] < double(convergenceCriteria)){
										calcFlagIn = 1;
									}

									if(jr > 0 && jth > 0 && jr < sampleSizeHalf && jth < sampleSizeHalf){
										if(ssfConvergenceIn[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] < double(convergenceCriteria)){
											calcFlagIn = 1;
										}
										if(ssfConvergenceIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] < double(convergenceCriteria)){
											calcFlagIn = 1;
										}
										if(ssfConvergenceIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + sampleSize - jth] < double(convergenceCriteria)){
											calcFlagIn = 1;
										}
									}else if(jr > 0 && jr < sampleSizeHalf){
										if(ssfConvergenceIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] < double(convergenceCriteria)){
											calcFlagIn = 1;
										}
									}else if(jth > 0 && jth < sampleSizeHalf){
										if(ssfConvergenceIn[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] < double(convergenceCriteria)){
											calcFlagIn = 1;
										}
									}
								}

                if(calcFlagIn == 1){
                  ssfModeTermIn = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, jth*polarSampleRate, In);

                  if( isnan(abs(ssfModeTermIn)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                    if(returnFlag < 2){
                      std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for In mode of F_"<<components[mu]<<". Set mode to zero. \n";
                      returnFlag = 2;
                    }
                  }else{
                    ssfIn[mu][l - abs(m)][jr*sampleSize + jth] += std::real(ssfModeTermIn*exp(-I*(phaseR + phaseTh)));
                    if(jr > 0 && jth > 0 && jr < sampleSizeHalf && jth < sampleSizeHalf){
                      ssfIn[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermIn*exp(-I*(phaseR - phaseTh)));
                      ssfIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermIn*exp(I*(phaseR - phaseTh)));
                      ssfIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + sampleSize - jth] += std::real(ssfModeTermIn*exp(I*(phaseR + phaseTh)));
                    }else if(jr > 0 && jr < sampleSizeHalf){
                      ssfIn[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermIn*exp(I*(phaseR - phaseTh)));
                    }else if(jth > 0 && jth < sampleSizeHalf){
                      ssfIn[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermIn*exp(-I*(phaseR - phaseTh)));
                    }
                  }
                }

								int calcFlagUp = 0;
								if(abs(teuk.getTeukolskyAmplitude(Up)) > 0.){
									if(ssfConvergenceUp[mu][l - abs(m)][jr*sampleSize + jth] < double(convergenceCriteria)){
										calcFlagUp = 1;
									}

									if(jr > 0 && jth > 0 && jr < sampleSizeHalf && jth < sampleSizeHalf){
										if(ssfConvergenceUp[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] < double(convergenceCriteria)){
											calcFlagUp = 1;
										}
										if(ssfConvergenceUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] < double(convergenceCriteria)){
											calcFlagUp = 1;
										}
										if(ssfConvergenceUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + sampleSize - jth] < double(convergenceCriteria)){
											calcFlagUp = 1;
										}
									}else if(jr > 0 && jr < sampleSizeHalf){
										if(ssfConvergenceUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] < double(convergenceCriteria)){
											calcFlagUp = 1;
										}
									}else if(jth > 0 && jth < sampleSizeHalf){
										if(ssfConvergenceUp[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] < double(convergenceCriteria)){
											calcFlagUp = 1;
										}
									}
								}

                if(calcFlagUp == 1){
                  ssfModeTermUp = self_force_amplitude(components[mu], teuk, l, m, jr*radialSampleRate, jth*polarSampleRate, Up);

                  if( isnan(abs(ssfModeTermUp)) || isnan(abs(phaseR)) || isnan(abs(phaseTh)) ){
                    if(returnFlag < 2){
                      std::cout << "(SSF) ERROR: (l,j,m,k,n) = ("<<l<<","<<j<<","<<m<<","<<k<<","<<n<<") returned NaN for Up mode of F_"<<components[mu]<<". Set mode to zero. \n";
                      returnFlag = 2;
                    }
                  }else{
                    ssfUp[mu][l - abs(m)][jr*sampleSize + jth] += std::real(ssfModeTermUp*exp(-I*(phaseR + phaseTh)));
                    if(jr > 0 && jth > 0 && jr < sampleSizeHalf && jth < sampleSizeHalf){
                      ssfUp[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermUp*exp(-I*(phaseR - phaseTh)));
                      ssfUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermUp*exp(I*(phaseR - phaseTh)));
                      ssfUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + sampleSize - jth] += std::real(ssfModeTermUp*exp(I*(phaseR + phaseTh)));
                    }else if(jr > 0 && jr < sampleSizeHalf){
                      ssfUp[mu][l - abs(m)][(sampleSize - jr)*sampleSize + jth] += std::real(ssfModeTermUp*exp(I*(phaseR - phaseTh)));
                    }else if(jth > 0 && jth < sampleSizeHalf){
                      ssfUp[mu][l - abs(m)][jr*sampleSize + sampleSize - jth] += std::real(ssfModeTermUp*exp(-I*(phaseR - phaseTh)));
                    }
                  }
                }

              }
            }
          }
        }
      }
    }
  }

  return returnFlag;
}

// int save_ssf_data_lm_ecceq(ComplexVector ssfData, int l, int m, GeodesicSource geo, const std::string &dir){
//   if(!boost::filesystem::exists(dir)){
//     std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
//     boost::filesystem::create_directory(dir);
//   }
//
//   char buff[7];
//   if(m < 0){
//     std::snprintf(buff, 7, "%d_n%d", l, -m);
//   }else{
//     std::snprintf(buff, 7, "%d_%d", l, m);
//   }
//   std::string filename = "ssf_data_";
// 	filename += buff;
// 	filename += ".txt";
// 	std::string filepath = dir + "/" + filename;
//
//   int Nr = geo.getRadialPosition().size() - 1;
//   int sampleSSF = sqrt(ssfData.size());
// 	if(!boost::filesystem::exists(filepath)){
// 		std::ofstream file;
// 		file.open(filepath);
// 		file << "r\tRe[Ft]\tIm[Ft]\n";
//
//     for(int jr = 0; jr <= sampleSSF/2; jr++){
//       file << std::scientific << std::setprecision(15);
//   		file << geo.getRadialPosition()[2*jr*Nr/sampleSSF] << "\t" << std::real(ssfData[jr*sampleSSF]) << "\t" << std::imag(ssfData[jr*sampleSSF]) << "\n";
//   	}
//
//   	for(int jr = sampleSSF/2 + 1; jr < sampleSSF; jr++){
//       file << std::scientific << std::setprecision(15);
//   		file << geo.getRadialPosition()[Nr - 2*(jr - sampleSSF/2)*Nr/sampleSSF] << "\t" << std::real(ssfData[jr*sampleSSF]) << "\t" << std::imag(ssfData[jr*sampleSSF]) << "\n";
//   	}
//
// 		file.close();
// 	}
//
//   return 0;
// }

int save_ssf_data(List components, SelfForceData ssfData, GeodesicSource geo, const std::string &dir){
	if(!boost::filesystem::exists(dir)){
		std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
		boost::filesystem::create_directory(dir);
	}
	std::string dirIn, dirUp;
	if(dir.back() == '/'){
		dirIn = dir + "in";
		dirUp = dir + "up";
	}else{
		dirIn = dir + "/in";
		dirUp = dir + "/up";
	}
  save_ssf_data(components, ssfData.in, 0, geo, dirIn, 0);
	save_ssf_data(components, ssfData.up, 0, geo, dirUp, 0);
	return 0;
}

int save_ssf_data(List components, SelfForceData ssfData, int m, GeodesicSource geo, const std::string &dir){
	if(!boost::filesystem::exists(dir)){
		std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
		boost::filesystem::create_directory(dir);
	}
	std::string dirIn, dirUp;
	if(dir.back() == '/'){
		dirIn = dir + "in";
		dirUp = dir + "up";
	}else{
		dirIn = dir + "/in";
		dirUp = dir + "/up";
	}
  save_ssf_data(components, ssfData.in, m, geo, dirIn, 1);
	save_ssf_data(components, ssfData.up, m, geo, dirUp, 1);
	return 0;
}

int save_ssf_data(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag){
  if(abs(geo.getInclination()) == 1. && geo.getEccentricity() == 0.){
    return save_ssf_data_circular(components, ssfData, m, geo, dir, lmFlag);
  }else if(abs(geo.getInclination()) == 1.){
    return save_ssf_data_equatorial(components, ssfData, m, geo, dir, lmFlag);
  }else if(geo.getEccentricity() == 0.){
    return save_ssf_data_spherical(components, ssfData, m, geo, dir, lmFlag);
  }
	return save_ssf_data_generic(components, ssfData, m, geo, dir, lmFlag);
}

int save_ssf_data_equatorial(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag){
  if(!boost::filesystem::exists(dir)){
    std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
    boost::filesystem::create_directory(dir);
  }

  int componentNum = ssfData.size();
  int lSize = ssfData[0].size();
  int Nr = geo.getRadialPosition().size() - 1;
  int sampleSSF = ssfData[0][0].size();

  for(int l = 0; l < lSize; l++){
    char buff[7];
		if(lmFlag){
			if(m < 0){
	      std::snprintf(buff, 7, "%d_n%d", l + abs(m), -m);
	    }else{
	      std::snprintf(buff, 7, "%d_%d", l + abs(m), m);
	    }
		}else{
			std::snprintf(buff, 7, "%d", l);
		}
    std::string filename = "ssf_data_";
  	filename += buff;
  	filename += ".txt";
		std::string filepath;
		if(dir.back() == '/'){
			filepath = dir + filename;
		}else{
			filepath = dir + "/" + filename;
		}

    std::ofstream file;
    file.open(filepath);
    file << "r\tth\tPhi\tFt\tFr\tFth\tFph\n";

    for(int jr = 0; jr <= sampleSSF/2; jr++){
      double Ft = 0.;
      double Fr = 0.;
      double Fth = 0.;
      double Fph = 0.;
			double phi = 0.;
      for(int mu = 0; mu < componentNum; mu++){
        if(components[mu] == 0){
          Ft = ssfData[mu][l][jr];
        }else if(components[mu] == 1){
          Fr = ssfData[mu][l][jr];
        }else if(components[mu] == 2){
          Fth = ssfData[mu][l][jr];
        }else if(components[mu] == 3){
          Fph = ssfData[mu][l][jr];
        }else if(components[mu] == 4){
					phi = ssfData[mu][l][jr];
				}
      }
      file << std::scientific << std::setprecision(15);
      file << geo.getRadialPosition(2*jr*Nr/sampleSSF) << "\t" << geo.getPolarPosition(0) << "\t" << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
    }

    for(int jr = sampleSSF/2 + 1; jr < sampleSSF; jr++){
      double Ft = 0.;
      double Fr = 0.;
      double Fth = 0.;
      double Fph = 0.;
			double phi = 0.;
      for(int mu = 0; mu < componentNum; mu++){
        if(components[mu] == 0){
          Ft = ssfData[mu][l][jr];
        }else if(components[mu] == 1){
          Fr = ssfData[mu][l][jr];
        }else if(components[mu] == 2){
          Fth = ssfData[mu][l][jr];
        }else if(components[mu] == 3){
          Fph = ssfData[mu][l][jr];
        }else if(components[mu] == 4){
					phi = ssfData[mu][l][jr];
				}
      }
      file << std::scientific << std::setprecision(15);
      file << geo.getRadialPosition(Nr - 2*(jr - sampleSSF/2)*Nr/sampleSSF) << "\t" << geo.getPolarPosition(0) << "\t" << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
    }

    file.close();
  }

  return 0;
}

int save_ssf_data_spherical(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag){
  if(!boost::filesystem::exists(dir)){
    std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
    boost::filesystem::create_directory(dir);
  }

  int componentNum = ssfData.size();
  int lSize = ssfData[0].size();
  int Nth = geo.getPolarPosition().size() - 1;
  int sampleSSF = ssfData[0][0].size();

  for(int l = 0; l < lSize; l++){
    char buff[7];
		if(lmFlag){
			if(m < 0){
	      std::snprintf(buff, 7, "%d_n%d", l + abs(m), -m);
	    }else{
	      std::snprintf(buff, 7, "%d_%d", l + abs(m), m);
	    }
		}else{
			std::snprintf(buff, 7, "%d", l);
		}
    std::string filename = "ssf_data_";
  	filename += buff;
  	filename += ".txt";
		std::string filepath;
		if(dir.back() == '/'){
			filepath = dir + filename;
		}else{
			filepath = dir + "/" + filename;
		}

    std::ofstream file;
    file.open(filepath);
    file << "r\tth\tPhi\tFt\tFr\tFth\tFph\n";

    for(int jth = 0; jth <= sampleSSF/2; jth++){
      double Ft = 0.;
      double Fr = 0.;
      double Fth = 0.;
      double Fph = 0.;
			double phi = 0.;
      for(int mu = 0; mu < componentNum; mu++){
        if(components[mu] == 0){
          Ft = ssfData[mu][l][jth];
        }else if(components[mu] == 1){
          Fr = ssfData[mu][l][jth];
        }else if(components[mu] == 2){
          Fth = ssfData[mu][l][jth];
        }else if(components[mu] == 3){
          Fph = ssfData[mu][l][jth];
        }else if(components[mu] == 4){
					phi = ssfData[mu][l][jth];
				}
      }
      file << std::scientific << std::setprecision(15);
      file << geo.getRadialPosition(0) << "\t" << geo.getPolarPosition(2*jth*Nth/sampleSSF) << "\t" << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
    }

    for(int jth = sampleSSF/2 + 1; jth < sampleSSF; jth++){
      double Ft = 0.;
      double Fr = 0.;
      double Fth = 0.;
      double Fph = 0.;
			double phi = 0.;
			for(int mu = 0; mu < componentNum; mu++){
        if(components[mu] == 0){
          Ft = ssfData[mu][l][jth];
        }else if(components[mu] == 1){
          Fr = ssfData[mu][l][jth];
        }else if(components[mu] == 2){
          Fth = ssfData[mu][l][jth];
        }else if(components[mu] == 3){
          Fph = ssfData[mu][l][jth];
        }else if(components[mu] == 4){
					phi = ssfData[mu][l][jth];
				}
      }
      file << std::scientific << std::setprecision(15);
      file << geo.getRadialPosition(0) << "\t" << geo.getPolarPosition(Nth - 2*(jth - sampleSSF/2)*Nth/sampleSSF) << "\t" << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
    }

    file.close();
  }

  return 0;
}

int save_ssf_data_circular(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag){
  if(!boost::filesystem::exists(dir)){
    std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
    boost::filesystem::create_directory(dir);
  }

  int componentNum = ssfData.size();
  int lSize = ssfData[0].size();

  for(int l = 0; l < lSize; l++){
    char buff[7];
		if(lmFlag){
			if(m < 0){
	      std::snprintf(buff, 7, "%d_n%d", l + abs(m), -m);
	    }else{
	      std::snprintf(buff, 7, "%d_%d", l + abs(m), m);
	    }
		}else{
			std::snprintf(buff, 7, "%d", l);
		}
    std::string filename = "ssf_data_";
  	filename += buff;
  	filename += ".txt";
		std::string filepath;
		if(dir.back() == '/'){
			filepath = dir + filename;
		}else{
			filepath = dir + "/" + filename;
		}

    std::ofstream file;
    file.open(filepath);
    file << "r\tth\tPhi\tFt\tFr\tFth\tFph\n";

		double Ft = 0.;
		double Fr = 0.;
		double Fth = 0.;
		double Fph = 0.;
		double phi = 0.;
		for(int mu = 0; mu < componentNum; mu++){
			if(components[mu] == 0){
				Ft = ssfData[mu][l][0];
			}else if(components[mu] == 1){
				Fr = ssfData[mu][l][0];
			}else if(components[mu] == 2){
				Fth = ssfData[mu][l][0];
			}else if(components[mu] == 3){
				Fph = ssfData[mu][l][0];
			}else if(components[mu] == 4){
				phi = ssfData[mu][l][0];
			}
		}
		file << std::scientific << std::setprecision(15);
		file << geo.getRadialPosition(0) << "\t" << geo.getPolarPosition(0) << "\t" << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";

    file.close();
  }

  return 0;
}

int save_ssf_data_generic(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag){
  if(!boost::filesystem::exists(dir)){
    std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
    boost::filesystem::create_directory(dir);
  }

  int componentNum = ssfData.size();
  int lSize = ssfData[0].size();
  int Nr = geo.getRadialPosition().size() - 1;
  int Nth = geo.getPolarPosition().size() - 1;
  int sampleSSF = sqrt(ssfData[0][0].size());

  for(int l = 0; l < lSize; l++){
    char buff[7];
		if(lmFlag){
			if(m < 0){
	      std::snprintf(buff, 7, "%d_n%d", l + abs(m), -m);
	    }else{
	      std::snprintf(buff, 7, "%d_%d", l + abs(m), m);
	    }
		}else{
			std::snprintf(buff, 7, "%d", l);
		}
    std::string filename = "ssf_data_";
  	filename += buff;
  	filename += ".txt";
		std::string filepath;
		if(dir.back() == '/'){
			filepath = dir + filename;
		}else{
			filepath = dir + "/" + filename;
		}

    std::ofstream file;
    file.open(filepath);
    file << "r\tth\tPhi\tFt\tFr\tFth\tFph\n";

    for(int jr = 0; jr <= sampleSSF/2; jr++){
      for(int jth = 0; jth <= sampleSSF/2; jth++){
        double Ft = 0.;
        double Fr = 0.;
        double Fth = 0.;
        double Fph = 0.;
        double phi = 0.;
        for(int mu = 0; mu < componentNum; mu++){
          if(components[mu] == 0){
            Ft = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 1){
            Fr = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 2){
            Fth = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 3){
            Fph = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 4){
            phi = ssfData[mu][l][jr*sampleSSF + jth];
          }
        }
        file << std::scientific << std::setprecision(15);
        file << geo.getRadialPosition()[2*jr*Nr/sampleSSF] << "\t" << geo.getPolarPosition()[2*jth*Nth/sampleSSF] << "\t"  << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
      }

      for(int jth = sampleSSF/2 + 1; jth < sampleSSF; jth++){
        double Ft = 0.;
        double Fr = 0.;
        double Fth = 0.;
        double Fph = 0.;
        double phi = 0.;
        for(int mu = 0; mu < componentNum; mu++){
          if(components[mu] == 0){
            Ft = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 1){
            Fr = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 2){
            Fth = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 3){
            Fph = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 4){
            phi = ssfData[mu][l][jr*sampleSSF + jth];
          }
        }
        file << std::scientific << std::setprecision(15);
        file << geo.getRadialPosition()[2*jr*Nr/sampleSSF] << "\t" << geo.getPolarPosition()[Nth - 2*(jth - sampleSSF/2)*Nth/sampleSSF]  << "\t"  << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
      }
      file << "\n";
    }

    for(int jr = sampleSSF/2 + 1; jr < sampleSSF; jr++){
      for(int jth = 0; jth <= sampleSSF/2; jth++){
        double Ft = 0.;
        double Fr = 0.;
        double Fth = 0.;
        double Fph = 0.;
        double phi = 0.;
        for(int mu = 0; mu < componentNum; mu++){
          if(components[mu] == 0){
            Ft = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 1){
            Fr = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 2){
            Fth = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 3){
            Fph = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 4){
            phi = ssfData[mu][l][jr*sampleSSF + jth];
          }
        }
        file << std::scientific << std::setprecision(15);
        file << geo.getRadialPosition()[Nr - 2*(jr - sampleSSF/2)*Nr/sampleSSF] << "\t" << geo.getPolarPosition()[2*jth*Nth/sampleSSF] << "\t" << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
      }

      for(int jth = sampleSSF/2 + 1; jth < sampleSSF; jth++){
        double Ft = 0.;
        double Fr = 0.;
        double Fth = 0.;
        double Fph = 0.;
        double phi = 0.;
        for(int mu = 0; mu < componentNum; mu++){
          if(components[mu] == 0){
            Ft = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 1){
            Fr = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 2){
            Fth = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 3){
            Fph = ssfData[mu][l][jr*sampleSSF + jth];
          }else if(components[mu] == 4){
            phi = ssfData[mu][l][jr*sampleSSF + jth];
          }
        }
        file << std::scientific << std::setprecision(15);
        file << geo.getRadialPosition()[Nr - 2*(jr - sampleSSF/2)*Nr/sampleSSF] << "\t" << geo.getPolarPosition()[Nth - 2*(jth - sampleSSF/2)*Nth/sampleSSF]  << "\t" << phi << "\t" << Ft << "\t" << Fr << "\t" << Fth << "\t" << Fph << "\n";
      }
      file << "\n";
    }

    file.close();
  }

  return 0;
}

int save_params(int lmax, int sampleNum, GeodesicSource& geo, const std::string &dir){
	std::string paramsfile;
	if(dir.back() == '/'){
		paramsfile = dir + "params.txt";
	}else{
		paramsfile = dir + "/params.txt";
	}
	std::ofstream file;
	file.open(paramsfile);
	file << "# Orbit data\n";
	file << "# a\tp\te\tx\tN\n";
	file << std::scientific << std::setprecision(15);
	file << geo.getBlackHoleSpin() << "\t" << geo.getSemiLatusRectum() << "\t" << geo.getEccentricity() << "\t" << geo.getInclination() << "\t" << geo.getRadialPosition().size() - 1 << "\n";
	file << "\n# Self-force data\n";
	file << "# lmax\tN\n";
	file << lmax << "\t" << sampleNum << "\n";
	file.close();
	return 0;
}

ParamsFileContents load_params_file(const std::string &dir){
	std::string paramsfile;
	if(dir.back() == '/'){
		paramsfile = dir + "params.txt";
	}else{
		paramsfile = dir + "/params.txt";
	}
	int lmaxTemp, orbitSampleNumTemp, sfSampleNumTemp;
	double aTemp, pTemp, eTemp, xTemp;
	ParamsFileContentsStruct params;

	std::istringstream lin;
	std::ifstream inFile(paramsfile);
	std::string line;
	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> aTemp >> pTemp >> eTemp >> xTemp >> orbitSampleNumTemp;
	params.a = aTemp;
	params.p = pTemp;
	params.e = eTemp;
	params.x = xTemp;
	params.orbitSampleNum = orbitSampleNumTemp;

	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty()){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> lmaxTemp >> sfSampleNumTemp;
	params.lmax = lmaxTemp;
	params.sfSampleNum = sfSampleNumTemp;
	inFile.close();

	return params;
}

SelfForceData load_ssf_data(const std::string &dir){
	ParamsFileContents params	= load_params_file(dir);
	int sfSampleSize = params.sfSampleNum;
	if(params.e == 0. && abs(params.x) == 1.){
		sfSampleSize = 1;
	}else if(params.e > 0. && abs(params.x) < 1.){
		sfSampleSize *= sfSampleSize;
	}
	int componentsNum = 5;
	RealTensor ssfDataUp = scalar_self_force_components_data_init(componentsNum, params.lmax, 0, sfSampleSize);
	RealTensor ssfDataIn = scalar_self_force_components_data_init(componentsNum, params.lmax, 0, sfSampleSize);
	double r, th;

	for(int l = 0; l < params.lmax; l++){
		char buff[7];
		std::snprintf(buff, 7, "%d", l);
		std::string filename = "ssf_data_";
		filename += buff;
		filename += ".txt";
		std::string filepathUp;
		std::string filepathIn;
		if(dir.back() == '/'){
			filepathUp = dir + "up/" + filename;
			filepathIn = dir + "in/" + filename;
		}else{
			filepathUp = dir + "/up/" + filename;
			filepathIn = dir + "/in/" + filename;
		}

		if(!boost::filesystem::exists(filepathUp)){
			std::cout << "Filepath " << filepathUp << " does not exist. \n";
		}
		std::istringstream lin;
		std::ifstream inFile(filepathUp);
		int i = 0;
		for(std::string line; std::getline(inFile, line); ) {
			lin.clear();
			lin.str(line);
			double Phi, Ft, Fr, Fth, Fph;
			if(!(line.front() == '#' || line.empty() || isalpha(line.front()))){
				lin >> r >> th >> Phi >> Ft >> Fr >> Fth >> Fph;
				ssfDataUp[0][l][i] = Ft;
				ssfDataUp[1][l][i] = Fr;
				ssfDataUp[2][l][i] = Fth;
				ssfDataUp[3][l][i] = Fph;
				ssfDataUp[4][l][i] = Phi;
				i++;
			}
		}
		inFile.close();
		if(i != sfSampleSize){
			std::cout << "(SSF) Error: Only " << i << " out of " << sfSampleSize << " points stored \n";
		}

		inFile.open(filepathIn);
		i = 0;
		for(std::string line; std::getline(inFile, line); ) {
			lin.clear();
			lin.str(line);
			double Phi, Ft, Fr, Fth, Fph;
			if(!(line.front() == '#' || line.empty() || isalpha(line.front()))){
				lin >> r >> th >> Phi >> Ft >> Fr >> Fth >> Fph;
				ssfDataIn[0][l][i] = Ft;
				ssfDataIn[1][l][i] = Fr;
				ssfDataIn[2][l][i] = Fth;
				ssfDataIn[3][l][i] = Fph;
				ssfDataIn[4][l][i] = Phi;
				i++;
			}
		}
		inFile.close();
		if(i != sfSampleSize){
			std::cout << "(SSF) Error: Only " << i << " out of " << sfSampleSize << " points stored \n";
		}
	}
	SelfForceData ssfData = {
		.in = ssfDataIn,
		.up = ssfDataUp
	};

	return ssfData;
}

SelfForceData load_ssf_data(int m, const std::string &dir){
	ParamsFileContents params	= load_params_file(dir);
	int sfSampleSize = params.sfSampleNum;
	int lsize = params.lmax + 1 - abs(m);
	if(params.e == 0. && abs(params.x) == 1.){
		sfSampleSize = 1;
	}else if(params.e > 0. && abs(params.x) < 1.){
		sfSampleSize *= sfSampleSize;
	}
	int componentsNum = 5;
	RealTensor ssfDataUp = scalar_self_force_components_data_init(componentsNum, params.lmax, m, sfSampleSize);
	RealTensor ssfDataIn = scalar_self_force_components_data_init(componentsNum, params.lmax, m, sfSampleSize);
	double r, th;

	for(int l = 0; l < lsize; l++){
		char buff[7];
		if(m < 0){
			std::snprintf(buff, 7, "%d_n%d", l + abs(m), -m);
		}else{
			std::snprintf(buff, 7, "%d_%d", l + abs(m), m);
		}
		std::string filename = "ssf_data_";
		filename += buff;
		filename += ".txt";
		std::string filepathUp;
		std::string filepathIn;
		if(dir.back() == '/'){
			filepathUp = dir + "up/" + filename;
			filepathIn = dir + "in/" + filename;
		}else{
			filepathUp = dir + "/up/" + filename;
			filepathIn = dir + "/in/" + filename;
		}

		if(!boost::filesystem::exists(filepathUp)){
			std::cout << "Filepath " << filepathUp << " does not exist. \n";
		}
		std::istringstream lin;
		std::ifstream inFile(filepathUp);
		int i = 0;
		for(std::string line; std::getline(inFile, line); ) {
			lin.clear();
			lin.str(line);
			double Phi, Ft, Fr, Fth, Fph;
			if(!(line.front() == '#' || line.empty() || isalpha(line.front()))){
				lin >> r >> th >> Phi >> Ft >> Fr >> Fth >> Fph;
				ssfDataUp[0][l][i] = Ft;
				ssfDataUp[1][l][i] = Fr;
				ssfDataUp[2][l][i] = Fth;
				ssfDataUp[3][l][i] = Fph;
				ssfDataUp[4][l][i] = Phi;
				i++;
			}
		}
		inFile.close();
		if(i != sfSampleSize){
			std::cout << "(SSF) Error: Only " << i << " out of " << sfSampleSize << " points stored \n";
		}

		inFile.open(filepathIn);
		i = 0;
		for(std::string line; std::getline(inFile, line); ) {
			lin.clear();
			lin.str(line);
			double Phi, Ft, Fr, Fth, Fph;
			if(!(line.front() == '#' || line.empty() || isalpha(line.front()))){
				lin >> r >> th >> Phi >> Ft >> Fr >> Fth >> Fph;
				ssfDataIn[0][l][i] = Ft;
				ssfDataIn[1][l][i] = Fr;
				ssfDataIn[2][l][i] = Fth;
				ssfDataIn[3][l][i] = Fph;
				ssfDataIn[4][l][i] = Phi;
				i++;
			}
		}
		inFile.close();
		if(i != sfSampleSize){
			std::cout << "(SSF) Error: Only " << i << " out of " << sfSampleSize << " points stored \n";
		}
	}
	SelfForceData ssfData = {
		.in = ssfDataIn,
		.up = ssfDataUp
	};

	return ssfData;
}

SelfForceData load_ssf_data(int l, int m, const std::string &dir){
	ParamsFileContents params	= load_params_file(dir);
	int sfSampleSize = params.sfSampleNum;
	if(params.e == 0. && abs(params.x) == 1.){
		sfSampleSize = 1;
	}else if(params.e > 0. && abs(params.x) < 1.){
		sfSampleSize *= sfSampleSize;
	}
	int componentsNum = 5;
	RealTensor ssfDataUp = scalar_self_force_components_data_init(componentsNum, m, m, sfSampleSize);
	RealTensor ssfDataIn = scalar_self_force_components_data_init(componentsNum, m, m, sfSampleSize);
	double r, th;

	char buff[7];
	if(m < 0){
		std::snprintf(buff, 7, "%d_n%d", l, -m);
	}else{
		std::snprintf(buff, 7, "%d_%d", l, m);
	}
	std::string filename = "ssf_data_";
	filename += buff;
	filename += ".txt";
	std::string filepathUp;
	std::string filepathIn;
	if(dir.back() == '/'){
		filepathUp = dir + "up/" + filename;
		filepathIn = dir + "in/" + filename;
	}else{
		filepathUp = dir + "/up/" + filename;
		filepathIn = dir + "/in/" + filename;
	}

	std::istringstream lin;
	std::ifstream inFile(filepathUp);
	int i = 0;
	for(std::string line; std::getline(inFile, line); ) {
		lin.clear();
		lin.str(line);
		double Phi, Ft, Fr, Fth, Fph;
		if(!(line.front() == '#' || line.empty() || isalpha(line.front()))){
			lin >> r >> th >> Phi >> Ft >> Fr >> Fth >> Fph;
			ssfDataUp[0][0][i] = Ft;
			ssfDataUp[1][0][i] = Fr;
			ssfDataUp[2][0][i] = Fth;
			ssfDataUp[3][0][i] = Fph;
			ssfDataUp[4][0][i] = Phi;
			i++;
		}
	}
	inFile.close();
	if(i != sfSampleSize){
		std::cout << "(SSF) Error: Only " << i << " out of " << sfSampleSize << " points stored \n";
	}
	inFile.open(filepathIn);
	i = 0;
	for(std::string line; std::getline(inFile, line); ) {
		lin.clear();
		lin.str(line);
		double Phi, Ft, Fr, Fth, Fph;
		if(!(line.front() == '#' || line.empty() || isalpha(line.front()))){
			lin >> r >> th >> Phi >> Ft >> Fr >> Fth >> Fph;
			ssfDataIn[0][0][i] = Ft;
			ssfDataIn[1][0][i] = Fr;
			ssfDataIn[2][0][i] = Fth;
			ssfDataIn[3][0][i] = Fph;
			ssfDataIn[4][0][i] = Phi;
			i++;
		}
	}
	inFile.close();
	if(i != sfSampleSize){
		std::cout << "(SSF) Error: Only " << i << " out of " << sfSampleSize << " points stored \n";
	}
	SelfForceData ssfData = {
		.in = ssfDataIn,
		.up = ssfDataUp
	};

	return ssfData;
}

SelfForceData load_and_construct_ssf_multipoles(const std::string &dir){
	std::cout << "(SSF) Loading m =  0 modes \n";
	SelfForceData ssf = load_ssf_data(0, dir);
	int lmax = ssf.in[0].size() - 1;
	// RealTensor ssfLIn = scalar_self_force_components_data_init(5, lmax, 0, ssf.in[0][0].size());
	// RealTensor ssfLUp = scalar_self_force_components_data_init(5, lmax, 0, ssf.in[0][0].size());
	// SelfForceData ssfL = {
	// 	.in = ssfLIn,
	// 	.up = ssfLUp
	// };
	SelfForceData ssfL = ssf;
	for(int m = 1; m <= lmax; m++){
		std::cout << "(SSF) Loading m =  "<<m<<" modes \n";
		ssf = load_ssf_data(m, dir);
		for(size_t mu = 0; mu < ssf.in.size(); mu++){
			for(size_t lm = 0; lm < ssf.in[mu].size(); lm++){
				for(size_t jj = 0; jj < ssf.in[mu][lm].size(); jj++){
					ssfL.in[mu][lm + abs(m)][jj] += ssf.in[mu][lm][jj];
					ssfL.up[mu][lm + abs(m)][jj] += ssf.up[mu][lm][jj];
				}
			}
		}
		std::cout << "(SSF) Loading m = "<<-m<<" modes \n";
		ssf = load_ssf_data(-m, dir);
		for(size_t mu = 0; mu < ssf.in.size(); mu++){
			for(size_t lm = 0; lm < ssf.in[mu].size(); lm++){
				for(size_t jj = 0; jj < ssf.in[mu][lm].size(); jj++){
					ssfL.in[mu][lm + abs(m)][jj] += ssf.in[mu][lm][jj];
					ssfL.up[mu][lm + abs(m)][jj] += ssf.up[mu][lm][jj];
				}
			}
		}
	}
	return ssfL;
}

int load_construct_and_save_ssf_multipoles(const std::string &dir){
	SelfForceData ssfL = load_and_construct_ssf_multipoles(dir);
	ParamsFileContents params = load_params_file(dir);
	char buff[500];
	if(dir.back() == '/'){
		sprintf(buff, "orbit/geo_a%.5f_p%.5f_e%.5f_x%.5f.txt", params.a, params.p, params.e, params.x);
	}else{
		sprintf(buff, "/orbit/geo_a%.5f_p%.5f_e%.5f_x%.5f.txt", params.a, params.p, params.e, params.x);
	}
	std::string geo_file = buff;
	GeodesicSource geo = load_geodesic(dir + geo_file);
	List components = {0, 1, 2, 3, 4};
	return save_ssf_data(components, ssfL, geo, dir);
}
