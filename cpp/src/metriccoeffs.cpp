// metriccoeffs.cpp

#include "metriccoeffs.hpp"

Complex metric_coefficient_ORG(int ai, int bi, int nt, int nr, int nz, int np, double a, double r, double z){
	if(ai > bi){
		return metric_coefficient_ORG(bi, ai, nt, nr, nz, np, a, r, z);
	}

	if(ai == 2 || bi == 2){
		return 0.;
	}else if(ai == 3 && bi == 4){
		return 0.;
	}

	if(ai == 1 && bi == 1){
		return metric_coefficient_ORG_11(nt, nr, nz, np, a, r, z);
	}else if(ai == 1 && bi == 3){
		return metric_coefficient_ORG_13(nt, nr, nz, np, a, r, z);
	}else if(ai == 1 && bi == 4){
		return std::conj(metric_coefficient_ORG_13(nt, nr, nz, np, a, r, z));
	}else if(ai == 3 && bi == 3){
		return metric_coefficient_ORG_33(nt, nr, nz, np, a, r, z);
	}else if(ai == 4 && bi == 4){
		return std::conj(metric_coefficient_ORG_33(nt, nr, nz, np, a, r, z));
	}else{
		std::cout << "(METRICCOEFFS) Error: Non-valid metric components (a, b) = ("<<ai<<", "<<bi<<")\n";
	}
	return 0.;
}

Complex metric_coefficient_IRG(int ai, int bi, int nt, int nr, int nz, int np, double a, double r, double z){
	if(ai > bi){
		return metric_coefficient_IRG(bi, ai, nt, nr, nz, np, a, r, z);
	}

	if(ai == 1 || bi == 1){
		return 0.;
	}else if(ai == 3 && bi == 4){
		return 0.;
	}

	if(ai == 2 && bi == 2){
		return metric_coefficient_IRG_22(nt, nr, nz, np, a, r, z);
	}else if(ai == 2 && bi == 3){
		return std::conj(metric_coefficient_IRG_24(nt, nr, nz, np, a, r, z));
	}else if(ai == 2 && bi == 4){
		return metric_coefficient_IRG_24(nt, nr, nz, np, a, r, z);
	}else if(ai == 3 && bi == 3){
		return std::conj(metric_coefficient_IRG_44(nt, nr, nz, np, a, r, z));
	}else if(ai == 4 && bi == 4){
		return metric_coefficient_IRG_44(nt, nr, nz, np, a, r, z);
	}else{
		std::cout << "(METRICCOEFFS) Error: Non-valid metric components (a, b) = ("<<ai<<", "<<bi<<")\n";
	}
	return 0.;
}

ComplexTensor metric_coefficients_ORG_11(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*(-1.*r[jr] + I*a*z[jz])*pow(1. - 1.*pow(z[jz],2),0.5);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.5*pow(r[jr] - I*a*z[jz],2);

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = (-1.*r[jr] + I*a*z[jz])*pow(a,2)*(-1. + pow(z[jz],2));

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*(-1.*r[jr] + I*a*z[jz])*(I*r[jr] + a*z[jz])*pow(1. - 1.*pow(z[jz],2),0.5);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.5*pow(a,2)*(-1. + pow(z[jz],2))*pow(r[jr] - I*a*z[jz],2);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_11_dz(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*(r[jr]*z[jz] + I*(a - 2.*a*pow(z[jz],2)))*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*(r[jr] - I*a*z[jz]);

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,2)*(-2.*r[jr]*z[jz] + I*a*(-1. + 3.*pow(z[jz],2)));

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*(I*z[jz]*pow(r[jr],2) + I*z[jz]*pow(a,2)*(2. - 3.*pow(z[jz],2)) + 2.*a*r[jr]*(-1. + 2.*pow(z[jz],2)))*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*(-1.*r[jr] + I*a*z[jz])*pow(a,2)*(-1.*r[jr]*z[jz] + I*a*(-1. + 2.*pow(z[jz],2)));
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_11_dz2(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*(r[jr] + I*a*z[jz]*(-3. + 2.*pow(z[jz],2)))*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,2);

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*(-1.*r[jr] + Complex(0.,3.)*a*z[jz])*pow(a,2);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*(a*r[jr]*z[jz]*(3. - 2.*pow(z[jz],2)) - 1.*a*r[jr]*z[jz]*(-3. + 2.*pow(z[jz],2)) + I*(pow(r[jr],2) + pow(a,2)*(2. - 9.*pow(z[jz],2) + 6.*pow(z[jz],4))))*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*pow(a,2)*(pow(r[jr],2) + I*a*(-6.*r[jr]*z[jz] + I*a*(-1. + 6.*pow(z[jz],2))));
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_13(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*(-1. + r[jr])*z[jz]*(r[jr] - I*a*z[jz])*pow(a,2)*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),0.5);

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = z[jz]*(r[jr] - I*a*z[jz])*pow(a,3)*(-1. + pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*(-1.*r[jr] + I*a*z[jz])*(pow(r[jr],3) + pow(a,2)*(-1.*r[jr] + 2.*(-1. + r[jr])*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.35355339059327373*a*(r[jr] - I*a*z[jz])*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = z[jz]*(-1.*r[jr] + I*a*z[jz])*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(-1. + pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.35355339059327373*(-1.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*(r[jr] - I*a*z[jz])*(-1. + pow(z[jz],2))*(a*z[jz]*(pow(a,2) + pow(r[jr],2)) + I*(pow(r[jr],3) + pow(a,2)*(-1.*r[jr] + 2.*(-1. + r[jr])*pow(z[jz],2))))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*(-1.*r[jr] + I*a*z[jz])*pow(a,2)*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.35355339059327373*(r[jr] - I*a*z[jz])*(pow(a,2) + pow(r[jr],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*a*(r[jr] - I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*a*(-1.*r[jr] + I*a*z[jz])*(pow(a,2) + pow(r[jr],2))*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_13_dz(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.8284271247461903*(-1. + r[jr])*pow(a,2)*(pow(r[jr],2)*(1. - 2.*pow(z[jz],2)) + Complex(0.,3.)*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*pow(a,3)*(pow(r[jr],2)*(1. - 2.*pow(z[jz],2)) + Complex(0.,3.)*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*a*(-4.*a*(-1. + r[jr])*z[jz]*pow(r[jr],2) + Complex(0.,3.)*r[jr]*(pow(r[jr],3) + pow(a,2)*(-1.*r[jr] + 2.*(-1. + r[jr])*pow(z[jz],2))) - 1.*z[jz]*(-1.*a*pow(r[jr],3) + pow(a,3)*(r[jr] + 2.*(-1. + r[jr])*pow(z[jz],2))))*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.35355339059327373*pow(a,2)*(-1.*a*z[jz]*(r[jr] - 1.*a*z[jz])*(r[jr] + a*z[jz]) + 2.*a*z[jz]*pow(r[jr],2) - Complex(0.,3.)*(pow(r[jr],3) + r[jr]*pow(a,2)*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(r[jr],2)*(1. - 2.*pow(z[jz],2)) + Complex(0.,3.)*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.35355339059327373*a*((-2. + r[jr])*r[jr] + pow(a,2))*(-1.*a*z[jz]*pow(r[jr],2) + Complex(0.,3.)*pow(r[jr],3) + Complex(0.,3.)*r[jr]*pow(a,2)*pow(z[jz],2) - 1.*pow(a,3)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*(I*z[jz]*pow(r[jr],5) + a*pow(r[jr],4)*(-1. + 3.*(-1. + pow(z[jz],2)) + 2.*pow(z[jz],2)) - I*z[jz]*pow(a,2)*pow(r[jr],2)*(-4. + r[jr]*(1. - 3.*pow(z[jz],2)) + 6.*pow(z[jz],2)) + pow(a,5)*pow(z[jz],4) + I*z[jz]*pow(a,4)*(2.*r[jr] + (-2.*(-1. + r[jr]) - 3.*r[jr])*pow(z[jz],2) + 4.*(-1. + r[jr])*pow(z[jz],4)) - 1.*r[jr]*pow(a,3)*(-2.*r[jr] + (-6. + 7.*r[jr])*pow(z[jz],2) - 1.*(-6. + 7.*r[jr])*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*pow(a,2)*(Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)*(-1. + 2.*pow(z[jz],2))) - 1.*z[jz]*pow(r[jr],2)*(pow(r[jr],2) + pow(a,2)*(-2. + 3.*pow(z[jz],2))))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.35355339059327373*a*(pow(a,2) + pow(r[jr],2))*(-1.*a*z[jz]*pow(r[jr],2) + Complex(0.,3.)*pow(r[jr],3) + Complex(0.,3.)*r[jr]*pow(a,2)*pow(z[jz],2) - 1.*pow(a,3)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = Complex(0.,-0.5)*a*((-2. + r[jr])*r[jr] + pow(a,2))*(Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)*(-1. + 2.*pow(z[jz],2))) - 1.*z[jz]*pow(r[jr],2)*(pow(r[jr],2) + pow(a,2)*(-2. + 3.*pow(z[jz],2))))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*a*(pow(a,2) + pow(r[jr],2))*(Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)*(-1. + 2.*pow(z[jz],2))) - 1.*z[jz]*pow(r[jr],2)*(pow(r[jr],2) + pow(a,2)*(-2. + 3.*pow(z[jz],2))))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_13_dz2(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.8284271247461903*(-1. + r[jr])*pow(a,2)*(-1.*r[jr]*z[jz]*pow(a,2)*(6. - 5.*pow(z[jz],2)) + z[jz]*pow(r[jr],3)*(-3. + 2.*pow(z[jz],2)) + I*a*pow(r[jr],2)*(-6. + 15.*pow(z[jz],2) - 10.*pow(z[jz],4)) - I*pow(a,3)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,3)*(-1.*r[jr]*z[jz]*pow(a,2)*(6. - 5.*pow(z[jz],2)) + z[jz]*pow(r[jr],3)*(-3. + 2.*pow(z[jz],2)) + I*a*pow(r[jr],2)*(-6. + 15.*pow(z[jz],2) - 10.*pow(z[jz],4)) - I*pow(a,3)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*r[jr]*pow(a,2)*(-5.*(a - 1.*r[jr])*r[jr]*(a + r[jr]) + Complex(0.,10.)*a*(-1. + r[jr])*r[jr]*z[jz] - 2.*(-1. + r[jr])*pow(r[jr],2) - I*a*z[jz]*(-1.*pow(a,2) + pow(r[jr],2)))*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*(-4.*r[jr] - 4.*I*a*z[jz])*pow(a,3)*pow(r[jr],2)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(-1.*r[jr]*z[jz]*pow(a,2)*(6. - 5.*pow(z[jz],2)) + z[jz]*pow(r[jr],3)*(-3. + 2.*pow(z[jz],2)) + I*a*pow(r[jr],2)*(-6. + 15.*pow(z[jz],2) - 10.*pow(z[jz],4)) - I*pow(a,3)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*(4.*r[jr] + 4.*I*a*z[jz])*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*pow(r[jr],2)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*(a*z[jz]*(pow(a,2) + pow(r[jr],2))*pow(r[jr],3)*(-3. + 2.*pow(z[jz],2)) + I*pow(r[jr],2)*(-1.*pow(r[jr],4) + pow(a,4)*(-6. + 15.*pow(z[jz],2) - 10.*pow(z[jz],4)) + r[jr]*pow(a,2)*(-4. - 1.*r[jr] - 3.*(-6. + r[jr])*pow(z[jz],2) + 2.*(-6. + r[jr])*pow(z[jz],4))) - 1.*a*r[jr]*z[jz]*(pow(r[jr],4)*(5. - 6.*pow(z[jz],2)) + pow(a,4)*(6. - 5.*pow(z[jz],2)) + r[jr]*pow(a,2)*(20. - 19.*r[jr] + (-42. + 43.*r[jr])*pow(z[jz],2) - 24.*(-1. + r[jr])*pow(z[jz],4))) - I*pow(a,2)*(pow(a,4)*pow(z[jz],4) + pow(r[jr],4)*(10. - 15.*pow(z[jz],2) + 6.*pow(z[jz],4)) + r[jr]*pow(a,2)*(-10.*r[jr] + 15.*r[jr]*pow(z[jz],2) + (-18. + 13.*r[jr])*pow(z[jz],4) - 16.*(-1. + r[jr])*pow(z[jz],6))) + z[jz]*pow(a,3)*(pow(r[jr],3)*(-2. + 3.*pow(z[jz],2)) + pow(a,2)*(2.*r[jr] - 3.*r[jr]*pow(z[jz],2) + 6.*(-1. + r[jr])*pow(z[jz],4) - 4.*(-1. + r[jr])*pow(z[jz],6))))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*pow(a,2)*(pow(r[jr],5) - 1.*pow(a,2)*pow(r[jr],3)*(-8. + 6.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],4)*(-5. + 6.*pow(z[jz],2)) - 1.*r[jr]*pow(a,4)*(-9. + 8.*pow(z[jz],2))*pow(z[jz],4) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 18.*pow(z[jz],2) + 12.*pow(z[jz],4)) - I*pow(a,5)*(-3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.7071067811865475*(4.*r[jr] + 4.*I*a*z[jz])*pow(a,2)*pow(r[jr],2)*(pow(a,2) + pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-4);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = Complex(0.,-0.5)*a*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(r[jr],5) - 1.*pow(a,2)*pow(r[jr],3)*(-8. + 6.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],4)*(-5. + 6.*pow(z[jz],2)) - 1.*r[jr]*pow(a,4)*(-9. + 8.*pow(z[jz],2))*pow(z[jz],4) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 18.*pow(z[jz],2) + 12.*pow(z[jz],4)) - I*pow(a,5)*(-3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*a*(pow(a,2) + pow(r[jr],2))*(pow(r[jr],5) - 1.*pow(a,2)*pow(r[jr],3)*(-8. + 6.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],4)*(-5. + 6.*pow(z[jz],2)) - 1.*r[jr]*pow(a,4)*(-9. + 8.*pow(z[jz],2))*pow(z[jz],4) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 18.*pow(z[jz],2) + 12.*pow(z[jz],4)) - I*pow(a,5)*(-3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_33(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = (r[jr] - I*a*z[jz])*(2.*r[jr] + I*a*(2. + 3.*(-2. + r[jr])*r[jr])*z[jz] + (-2. + r[jr])*pow(a,2) + I*z[jz]*pow(a,3) - 1.*pow(r[jr],3))*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*a*(-1.*r[jr] + I*a*z[jz])*(r[jr] + Complex(0.,3.)*a*(-1. + r[jr])*z[jz] + pow(a,2) - 2.*pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.25*pow(a,2)*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*(r[jr] - I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*((2. - 3.*r[jr])*r[jr] + 4.*I*a*(-1. + r[jr])*z[jz] + pow(a,2))*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*(a*(-2. + r[jr])*r[jr] + pow(a,3))*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.25*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),2);

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*(r[jr] - I*a*z[jz])*(r[jr]*(-1. + 2.*r[jr])*pow(a,2) + I*(3. - 4.*r[jr])*z[jz]*pow(a,3) - 1.*pow(a,4) + I*a*(5. - 4.*r[jr])*z[jz]*pow(r[jr],2) + 3.*(-1. + r[jr])*pow(r[jr],3))*pow(r[jr] + I*a*z[jz],-2);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.5*a*(pow(a,2) + pow(r[jr],2))*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2))*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.25*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow(pow(a,2) + pow(r[jr],2),2);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_33_dz(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = Complex(0.,-2.)*a*(I*a*r[jr]*(2. + r[jr]*(-9. + 5.*r[jr]))*z[jz] + (-3. + r[jr])*r[jr]*pow(a,2) + I*(1. + r[jr])*z[jz]*pow(a,3) + (2. - 3.*(-1. + r[jr])*r[jr])*pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*pow(a,2)*((-I)*a*z[jz]*((10. - 11.*r[jr])*r[jr] + pow(a,2)) + 3.*r[jr]*((2. - 3.*r[jr])*r[jr] + pow(a,2)))*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*r[jr]*(r[jr] - I*a*z[jz])*pow(a,3)*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*a*((-2. + r[jr])*r[jr] + pow(a,2))*(r[jr]*(r[jr]*(-10. + 13.*r[jr]) - 3.*pow(a,2)) + I*a*z[jz]*((14. - 15.*r[jr])*r[jr] + pow(a,2)))*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*I*r[jr]*(-1.*r[jr] + I*a*z[jz])*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = (-I)*a*r[jr]*(-1.*r[jr] + I*a*z[jz])*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),2);

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*a*(r[jr]*(2.*(3. - 5.*r[jr])*r[jr]*pow(a,2) + 3.*pow(a,4) + (14. - 13.*r[jr])*pow(r[jr],3)) - I*a*z[jz]*(2.*(5. - 7.*r[jr])*r[jr]*pow(a,2) + pow(a,4) + 3.*(6. - 5.*r[jr])*pow(r[jr],3)))*pow(r[jr] + I*a*z[jz],-3);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = Complex(0.,-2.)*r[jr]*(-1.*r[jr] + I*a*z[jz])*pow(a,2)*(pow(a,2) + pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-3);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*I*a*r[jr]*(-1.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-3);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = (-I)*a*r[jr]*(-1.*r[jr] + I*a*z[jz])*pow(r[jr] + I*a*z[jz],-3)*pow(pow(a,2) + pow(r[jr],2),2);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_ORG_33_dz2(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -4.*pow(a,2)*(I*a*r[jr]*(2. + r[jr]*(-9. + 5.*r[jr]))*z[jz] + (-5. + r[jr])*r[jr]*pow(a,2) + I*(1. + r[jr])*z[jz]*pow(a,3) + (2. + (9. - 7.*r[jr])*r[jr])*pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*(r[jr]*(r[jr]*(-14. + 19.*r[jr]) - 5.*pow(a,2)) + I*a*z[jz]*((10. - 11.*r[jr])*r[jr] + pow(a,2)))*pow(a,3)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -2.*r[jr]*(-2.*r[jr] + I*a*z[jz])*pow(a,4)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(r[jr]*(r[jr]*(-22. + 27.*r[jr]) - 5.*pow(a,2)) + I*a*z[jz]*((14. - 15.*r[jr])*r[jr] + pow(a,2)))*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 4.*r[jr]*(-2.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*pow(a,3)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -2.*r[jr]*(-2.*r[jr] + I*a*z[jz])*pow(a,2)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),2);

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*pow(a,2)*(I*a*z[jz]*(2.*(5. - 7.*r[jr])*r[jr]*pow(a,2) + pow(a,4) + 3.*(6. - 5.*r[jr])*pow(r[jr],3)) + r[jr]*(2.*r[jr]*(-7. + 11.*r[jr])*pow(a,2) - 5.*pow(a,4) + 3.*(-10. + 9.*r[jr])*pow(r[jr],3)))*pow(r[jr] + I*a*z[jz],-4);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -4.*r[jr]*(-2.*r[jr] + I*a*z[jz])*pow(a,3)*(pow(a,2) + pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-4);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 4.*r[jr]*(-2.*r[jr] + I*a*z[jz])*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-4);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -2.*r[jr]*(-2.*r[jr] + I*a*z[jz])*pow(a,2)*pow(r[jr] + I*a*z[jz],-4)*pow(pow(a,2) + pow(r[jr],2),2);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_22(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),0.5);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.5*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,2)*(-1. + pow(z[jz],2))*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*(I*r[jr] + a*z[jz])*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),0.5);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.5*pow(a,2)*(-1. + pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_22_dz(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*(z[jz]*pow(r[jr],2) - I*a*r[jr]*(-1. + pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(-3. + 2.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,2)*(2.*z[jz]*pow(r[jr],2) - 1.*z[jz]*pow(a,2)*(-3. + pow(z[jz],2)) - I*a*r[jr]*(-1. + pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*((-I)*z[jz]*pow(r[jr],3) + a*pow(r[jr],2)*(2. - 3.*pow(z[jz],2)) + pow(a,3)*(-2. + pow(z[jz],2))*pow(z[jz],2) + I*r[jr]*z[jz]*pow(a,2)*(-4. + 3.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*(I*a + r[jr]*z[jz])*pow(a,2)*pow(r[jr] + I*a*z[jz],-3);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_22_dz2(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = (-I)*a*(pow(r[jr],4) + 2.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(a,4)*pow(z[jz],2)*(-12. + 19.*pow(z[jz],2) - 6.*pow(z[jz],4)) - 2.*I*r[jr]*z[jz]*pow(a,3)*(4. - 7.*pow(z[jz],2) + 3.*pow(z[jz],4)) + 2.*pow(a,2)*pow(r[jr],2)*(2. - 6.*pow(z[jz],2) + 5.*pow(z[jz],4)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 3.*pow(a,2)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*pow(a,2)*(Complex(0.,-2.)*a*z[jz]*pow(r[jr],3) + pow(r[jr],4) + 2.*I*r[jr]*z[jz]*pow(a,3)*(-2. + pow(z[jz],2)) + pow(a,4)*(-6. + pow(z[jz],2))*pow(z[jz],2) - 2.*pow(a,2)*pow(r[jr],2)*(-1. + 3.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*a*(a*z[jz]*pow(r[jr],4)*(-3. + 2.*pow(z[jz],2)) + I*r[jr]*pow(a,4)*pow(z[jz],2)*(12. - 19.*pow(z[jz],2) + 6.*pow(z[jz],4)) + z[jz]*pow(a,3)*(pow(a,2)*pow(z[jz],2)*(6. - 9.*pow(z[jz],2) + 2.*pow(z[jz],4)) - 2.*pow(r[jr],2)*(4. - 7.*pow(z[jz],2) + 3.*pow(z[jz],4))) - I*pow(r[jr],3)*(pow(r[jr],2) + pow(a,2)*(2. - 6.*pow(z[jz],2) + 4.*pow(z[jz],4))) - 2.*I*r[jr]*pow(a,2)*(pow(a,2)*pow(z[jz],2)*(-3. + 5.*pow(z[jz],2) - 2.*pow(z[jz],4)) + pow(r[jr],2)*(2. - 6.*pow(z[jz],2) + 5.*pow(z[jz],4))) - 2.*a*z[jz]*pow(r[jr],2)*(-1.*pow(r[jr],2)*(-1. + pow(z[jz],2)) + pow(a,2)*(5. - 10.*pow(z[jz],2) + 6.*pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = pow(a,2)*(I*a*(Complex(0.,3.)*a + 2.*r[jr]*z[jz]) - 1.*pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-4);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_24(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = z[jz]*pow(a,3)*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),0.5);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*r[jr]*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.7071067811865475*a*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*z[jz]*pow(a,2)*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),0.5);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2);

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = a*((-I)*r[jr]*((-2. + r[jr])*r[jr] + pow(a,2)) + a*z[jz]*(pow(a,2) + pow(r[jr],2)))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),0.5);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.5*I*pow(a,2)*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),0.5);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -0.7071067811865475*(pow(a,2) + pow(r[jr],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = Complex(0.,-0.5)*a*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),0.5);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = (-I)*a*(pow(a,2) + pow(r[jr],2))*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_24_dz(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*pow(a,3)*((-I)*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*(-2. + pow(z[jz],2))*pow(z[jz],2) + pow(r[jr],2)*(-1. + 2.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = Complex(0.,1.4142135623730951)*a*r[jr]*(r[jr] - Complex(0.,3.)*a*z[jz])*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*pow(a,2)*(-1.*a*z[jz]*pow(r[jr],2) - I*pow(r[jr],3) - I*r[jr]*pow(a,2)*pow(z[jz],2) - 1.*pow(a,3)*pow(z[jz],3))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*pow(a,2)*((-I)*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*(-2. + pow(z[jz],2))*pow(z[jz],2) + pow(r[jr],2)*(-1. + 2.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*a*(-1.*a*z[jz]*pow(r[jr],2) - I*pow(r[jr],3) - I*r[jr]*pow(a,2)*pow(z[jz],2) - 1.*pow(a,3)*pow(z[jz],3))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3);

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*a*((-I)*(-2. + r[jr])*z[jz]*pow(r[jr],4) - I*r[jr]*z[jz]*pow(a,4)*(2. - 1.*pow(z[jz],2)) - 1.*pow(a,5)*(-2. + pow(z[jz],2))*pow(z[jz],2) - I*z[jz]*pow(a,2)*pow(r[jr],2)*(r[jr]*pow(z[jz],2) - 1.*(-2. + r[jr])*(-3. + 2.*pow(z[jz],2))) + a*pow(r[jr],3)*(-1.*(-2. + r[jr])*(-1. + pow(z[jz],2)) + r[jr]*(-1. + 2.*pow(z[jz],2))) + pow(a,3)*pow(r[jr],2)*(3.*pow(z[jz],2) - 1.*pow(z[jz],4)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-0.5);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*pow(a,2)*((-I)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) + z[jz]*pow(a,2)*(pow(r[jr],2)*(3. - 2.*pow(z[jz],2)) + pow(a,2)*pow(z[jz],2)) + z[jz]*pow(r[jr],2)*(pow(r[jr],2) + pow(a,2)*(-2. + 3.*pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.7071067811865475*a*(pow(a,2) + pow(r[jr],2))*(-1.*a*z[jz]*pow(r[jr],2) - I*pow(r[jr],3) - I*r[jr]*pow(a,2)*pow(z[jz],2) - 1.*pow(a,3)*pow(z[jz],3))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*((-I)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) + z[jz]*pow(a,2)*(pow(r[jr],2)*(3. - 2.*pow(z[jz],2)) + pow(a,2)*pow(z[jz],2)) + z[jz]*pow(r[jr],2)*(pow(r[jr],2) + pow(a,2)*(-2. + 3.*pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = (-I)*a*(pow(a,2) + pow(r[jr],2))*(I*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2)*(3. - 2.*pow(z[jz],2)) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(r[jr],2)*(pow(r[jr],2) + pow(a,2)*(-2. + 3.*pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_24_dz2(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*pow(a,3)*(z[jz]*pow(r[jr],4)*(-3. + 2.*pow(z[jz],2)) + pow(a,4)*pow(z[jz],3)*(6. - 9.*pow(z[jz],2) + 2.*pow(z[jz],4)) + 2.*I*r[jr]*pow(a,3)*pow(z[jz],2)*(3. - 5.*pow(z[jz],2) + 2.*pow(z[jz],4)) - 2.*I*a*pow(r[jr],3)*(1. - 3.*pow(z[jz],2) + 2.*pow(z[jz],4)) - 2.*z[jz]*pow(a,2)*pow(r[jr],2)*(5. - 10.*pow(z[jz],2) + 6.*pow(z[jz],4)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 5.656854249492381*r[jr]*pow(a,2)*(Complex(0.,-2.)*a*r[jr]*z[jz] + pow(r[jr],2) - 3.*pow(a,2)*pow(z[jz],2))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.4142135623730951*pow(a,3)*(2.*I*a*z[jz]*pow(r[jr],3) - 1.*pow(r[jr],4) + 2.*I*r[jr]*pow(a,3)*pow(z[jz],3) + pow(a,4)*pow(z[jz],4))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*pow(a,2)*(z[jz]*pow(r[jr],4)*(-3. + 2.*pow(z[jz],2)) + pow(a,4)*pow(z[jz],3)*(6. - 9.*pow(z[jz],2) + 2.*pow(z[jz],4)) + 2.*I*r[jr]*pow(a,3)*pow(z[jz],2)*(3. - 5.*pow(z[jz],2) + 2.*pow(z[jz],4)) - 2.*I*a*pow(r[jr],3)*(1. - 3.*pow(z[jz],2) + 2.*pow(z[jz],4)) - 2.*z[jz]*pow(a,2)*pow(r[jr],2)*(5. - 10.*pow(z[jz],2) + 6.*pow(z[jz],4)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 1.4142135623730951*pow(a,2)*(2.*I*a*z[jz]*pow(r[jr],3) - 1.*pow(r[jr],4) + 2.*I*r[jr]*pow(a,3)*pow(z[jz],3) + pow(a,4)*pow(z[jz],4))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4);

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.4142135623730951*a*(-1.*a*z[jz]*(pow(a,2) + pow(r[jr],2))*pow(r[jr],4)*(-3. + 2.*pow(z[jz],2)) + I*r[jr]*((-2. + r[jr])*r[jr] + pow(a,2))*pow(a,4)*pow(z[jz],2)*(12. - 19.*pow(z[jz],2) + 6.*pow(z[jz],4)) + I*pow(r[jr],3)*(-1.*(-2. + r[jr])*pow(r[jr],3) + pow(a,2)*pow(r[jr],2)*(1. - 6.*pow(z[jz],2) + 4.*pow(z[jz],4)) + pow(a,4)*(2. - 6.*pow(z[jz],2) + 4.*pow(z[jz],4))) + 2.*a*z[jz]*pow(r[jr],2)*((-2. + r[jr])*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(a,4)*(5. - 10.*pow(z[jz],2) + 6.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],2)*(4. - 9.*pow(z[jz],2) + 6.*pow(z[jz],4))) + z[jz]*pow(a,3)*(pow(a,4)*pow(z[jz],2)*(-6. + 9.*pow(z[jz],2) - 2.*pow(z[jz],4)) - 2.*(-2. + r[jr])*pow(r[jr],3)*(4. - 7.*pow(z[jz],2) + 3.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],2)*(-8. + 8.*pow(z[jz],2) + 3.*pow(z[jz],4) - 2.*pow(z[jz],6))) - 2.*I*pow(a,2)*(r[jr]*pow(a,4)*pow(z[jz],2)*(3. - 5.*pow(z[jz],2) + 2.*pow(z[jz],4)) + (-2. + r[jr])*pow(r[jr],4)*(2. - 6.*pow(z[jz],2) + 5.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],3)*(2. - 3.*pow(z[jz],2) + 2.*pow(z[jz],6))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-1.5);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*pow(a,2)*(pow(r[jr],6) + 4.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(-1. + pow(z[jz],2)) + 2.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) + pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(19.*pow(z[jz],2) + 10.*(-2. + pow(z[jz],2))*pow(z[jz],2) - 6.*pow(z[jz],4)) + pow(a,6)*(-2. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 3.*pow(z[jz],2) + pow(z[jz],4)) - 1.*pow(a,2)*pow(r[jr],4)*(2. - 9.*pow(z[jz],2) + 6.*pow(z[jz],4) - 2.*(2. - 6.*pow(z[jz],2) + 5.*pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.4142135623730951*pow(a,2)*(pow(a,2) + pow(r[jr],2))*(2.*I*a*z[jz]*pow(r[jr],3) - 1.*pow(r[jr],4) + 2.*I*r[jr]*pow(a,3)*pow(z[jz],3) + pow(a,4)*pow(z[jz],4))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*(pow(r[jr],6) + 4.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(-1. + pow(z[jz],2)) + 2.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) + pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(19.*pow(z[jz],2) + 10.*(-2. + pow(z[jz],2))*pow(z[jz],2) - 6.*pow(z[jz],4)) + pow(a,6)*(-2. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 3.*pow(z[jz],2) + pow(z[jz],4)) - 1.*pow(a,2)*pow(r[jr],4)*(2. - 9.*pow(z[jz],2) + 6.*pow(z[jz],4) - 2.*(2. - 6.*pow(z[jz],2) + 5.*pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = I*a*(pow(a,2) + pow(r[jr],2))*(pow(r[jr],6) + 4.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(-1. + pow(z[jz],2)) + 2.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) + pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(19.*pow(z[jz],2) + 10.*(-2. + pow(z[jz],2))*pow(z[jz],2) - 6.*pow(z[jz],4)) + pow(a,6)*(-2. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 3.*pow(z[jz],2) + pow(z[jz],4)) - 1.*pow(a,2)*pow(r[jr],4)*(2. - 9.*pow(z[jz],2) + 6.*pow(z[jz],4) - 2.*(2. - 6.*pow(z[jz],2) + 5.*pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_44(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -2.*a*(r[jr]*(-3. + 2.*r[jr]) - I*a*(-1. + r[jr])*z[jz] + pow(a,2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*pow(a,2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2);

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*pow(r[jr] - I*a*z[jz],-1);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -2.*a*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*(-1.*pow(a,2) + pow(r[jr],2) - 1.*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2))*pow(-1.*r[jr] + I*a*z[jz],-1))*pow((-2. + r[jr])*r[jr] + pow(a,2),-2);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -2.*a*(pow(a,2) + pow(r[jr],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-2);

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -2. - 4.*r[jr]*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -1.*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(pow(a,2) + pow(r[jr],2),2);
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_44_dz(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*I*pow(a,2)*pow(r[jr] - I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*I*a*pow(r[jr] - I*a*z[jz],-2);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 2.*I*a*(pow(a,2) + pow(r[jr],2))*pow(r[jr] - I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;
		}
	}

	return coeffs;
}

ComplexTensor metric_coefficients_IRG_44_dz2(double a, Vector r, Vector z){
	int componentNum = 3*3*3*3;
	int ai, bi, ci, di, pos;
	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

	for(size_t jr = 0; jr < r.size(); jr++){
		for(size_t jz = 0; jz < z.size(); jz++){
		ai = 0, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 4.*pow(a,3)*pow(-1.*r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 0, bi = 0, ci = 0, di = 2;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 1, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 0, ci = 2, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = -4.*pow(a,2)*pow(r[jr] - I*a*z[jz],-3);

		ai = 0, bi = 1, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 1, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 0, bi = 2, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 4.*pow(a,2)*(pow(a,2) + pow(r[jr],2))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1);

		ai = 1, bi = 0, ci = 0, di = 1;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 0, ci = 1, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 1, bi = 1, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;

		ai = 2, bi = 0, ci = 0, di = 0;
		pos = 3*(3*((3*ai) + bi) + ci) + di;
		coeffs[pos][jr][jz] = 0.;
		}
	}

	return coeffs;
}

Complex metric_coefficient_IRG_22(int nt, int nr, int nz, int np, double a, double r, double z){
	Complex hab = 0.;

	if(nt == 0&& nr == 0&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 2){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 0){
		hab = I*a*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow(1. - 1.*pow(z,2),0.5);
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 2&& np == 0){
		hab = -0.5*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 2&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 0){
		hab = pow(a,2)*(-1. + pow(z,2))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2);
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 1&& nr == 0&& nz == 1&& np == 0){
		hab = a*(I*r + a*z)*pow(r - I*a*z,-1)*pow(r + I*a*z,-2)*pow(1. - 1.*pow(z,2),0.5);
	}else if(nt == 1&& nr == 1&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 2&& nr == 0&& nz == 0&& np == 0){
		hab = -0.5*pow(a,2)*(-1. + pow(z,2))*pow(r + I*a*z,-2);
	}

	return hab;
}

Complex metric_coefficient_IRG_24(int nt, int nr, int nz, int np, double a, double r, double z){
	Complex hab = 0.;

	if(nt == 0&& nr == 0&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 1){
		hab = z*pow(a,3)*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),0.5);
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 2){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 0){
		hab = 1.4142135623730951*r*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 1){
		hab = -0.7071067811865475*a*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1);
	}else if(nt == 0&& nr == 0&& nz == 2&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 0){
		hab = -1.*z*pow(a,2)*pow(r - I*a*z,-1)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),0.5);
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 1&& np == 0){
		hab = 0.7071067811865475*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 2&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 0){
		hab = a*((-I)*r*((-2. + r)*r + pow(a,2)) + a*z*(pow(a,2) + pow(r,2)))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),0.5);
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 1){
		hab = 0.5*I*pow(a,2)*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),0.5);
	}else if(nt == 1&& nr == 0&& nz == 1&& np == 0){
		hab = -0.7071067811865475*(pow(a,2) + pow(r,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1);
	}else if(nt == 1&& nr == 1&& nz == 0&& np == 0){
		hab = Complex(0.,-0.5)*a*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),0.5);
	}else if(nt == 2&& nr == 0&& nz == 0&& np == 0){
		hab = (-I)*a*(pow(a,2) + pow(r,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),-0.5);
	}

	return hab;
}

Complex metric_coefficient_IRG_44(int nt, int nr, int nz, int np, double a, double r, double z){
	Complex hab = 0.;

	if(nt == 0&& nr == 0&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 1){
		hab = -2.*a*(r*(-3. + 2.*r) - I*a*(-1. + r)*z + pow(a,2))*pow(-1.*r + I*a*z,-1)*pow((-2. + r)*r + pow(a,2),-2);
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 2){
		hab = -1.*pow(a,2)*pow((-2. + r)*r + pow(a,2),-2);
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 2&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 0){
		hab = 2.*pow(r - I*a*z,-1);
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 1){
		hab = -2.*a*pow((-2. + r)*r + pow(a,2),-1);
	}else if(nt == 0&& nr == 1&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 2&& nz == 0&& np == 0){
		hab = -1.;
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 0){
		hab = 2.*(-1.*pow(a,2) + pow(r,2) - 1.*((-2. + r)*r + pow(a,2))*(pow(a,2) + pow(r,2))*pow(-1.*r + I*a*z,-1))*pow((-2. + r)*r + pow(a,2),-2);
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 1){
		hab = -2.*a*(pow(a,2) + pow(r,2))*pow((-2. + r)*r + pow(a,2),-2);
	}else if(nt == 1&& nr == 0&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 1&& nr == 1&& nz == 0&& np == 0){
		hab = -2. - 4.*r*pow((-2. + r)*r + pow(a,2),-1);
	}else if(nt == 2&& nr == 0&& nz == 0&& np == 0){
		hab = -1.*pow((-2. + r)*r + pow(a,2),-2)*pow(pow(a,2) + pow(r,2),2);
	}

	return hab;
}

Complex metric_coefficient_ORG_11(int nt, int nr, int nz, int np, double a, double r, double z){
	Complex hab = 0.;

	if(nt == 0&& nr == 0&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 2){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 0){
		hab = I*a*(-1.*r + I*a*z)*pow(1. - 1.*pow(z,2),0.5);
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 2&& np == 0){
		hab = -0.5*pow(r - I*a*z,2);
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 2&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 0){
		hab = (-1.*r + I*a*z)*pow(a,2)*(-1. + pow(z,2));
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 1&& nr == 0&& nz == 1&& np == 0){
		hab = a*(-1.*r + I*a*z)*(I*r + a*z)*pow(1. - 1.*pow(z,2),0.5);
	}else if(nt == 1&& nr == 1&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 2&& nr == 0&& nz == 0&& np == 0){
		hab = -0.5*pow(a,2)*(-1. + pow(z,2))*pow(r - I*a*z,2);
	}

	return hab;
}

Complex metric_coefficient_ORG_13(int nt, int nr, int nz, int np, double a, double r, double z){
	Complex hab = 0.;

	if(nt == 0&& nr == 0&& nz == 0&& np == 0){
		hab = 2.*(-1. + r)*z*(r - I*a*z)*pow(a,2)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),0.5);
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 1){
		hab = z*(r - I*a*z)*pow(a,3)*(-1. + pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 2){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 0){
		hab = 0.7071067811865475*(-1.*r + I*a*z)*(pow(r,3) + pow(a,2)*(-1.*r + 2.*(-1. + r)*pow(z,2)))*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 1){
		hab = 0.35355339059327373*a*(r - I*a*z)*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 0&& nz == 2&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 0){
		hab = z*(-1.*r + I*a*z)*pow(a,2)*((-2. + r)*r + pow(a,2))*(-1. + pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 1&& np == 0){
		hab = 0.35355339059327373*(-1.*r + I*a*z)*((-2. + r)*r + pow(a,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 2&& nz == 0&& np == 0){
		hab = 0.;
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 0){
		hab = a*(r - I*a*z)*(-1. + pow(z,2))*(a*z*(pow(a,2) + pow(r,2)) + I*(pow(r,3) + pow(a,2)*(-1.*r + 2.*(-1. + r)*pow(z,2))))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 1){
		hab = 0.5*I*(-1.*r + I*a*z)*pow(a,2)*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
	}else if(nt == 1&& nr == 0&& nz == 1&& np == 0){
		hab = 0.35355339059327373*(r - I*a*z)*(pow(a,2) + pow(r,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2);
	}else if(nt == 1&& nr == 1&& nz == 0&& np == 0){
		hab = 0.5*I*a*(r - I*a*z)*((-2. + r)*r + pow(a,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
	}else if(nt == 2&& nr == 0&& nz == 0&& np == 0){
		hab = 0.5*I*a*(-1.*r + I*a*z)*(pow(a,2) + pow(r,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
	}

	return hab;
}

Complex metric_coefficient_ORG_33(int nt, int nr, int nz, int np, double a, double r, double z){
	Complex hab = 0.;

	if(nt == 0&& nr == 0&& nz == 0&& np == 0){
		hab = (r - I*a*z)*(2.*r + I*a*(2. + 3.*(-2. + r)*r)*z + (-2. + r)*pow(a,2) + I*z*pow(a,3) - 1.*pow(r,3))*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 1){
		hab = 0.5*a*(-1.*r + I*a*z)*(r + Complex(0.,3.)*a*(-1. + r)*z + pow(a,2) - 2.*pow(r,2))*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 0&& nz == 0&& np == 2){
		hab = -0.25*pow(a,2)*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 1&& np == 1){
		hab = 0.;
	}else if(nt == 0&& nr == 0&& nz == 2&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 0){
		hab = 0.5*(r - I*a*z)*((-2. + r)*r + pow(a,2))*((2. - 3.*r)*r + 4.*I*a*(-1. + r)*z + pow(a,2))*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 1&& nz == 0&& np == 1){
		hab = 0.5*(a*(-2. + r)*r + pow(a,3))*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
	}else if(nt == 0&& nr == 1&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 0&& nr == 2&& nz == 0&& np == 0){
		hab = -0.25*pow(r - I*a*z,2)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),2);
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 0){
		hab = 0.5*(r - I*a*z)*(r*(-1. + 2.*r)*pow(a,2) + I*(3. - 4.*r)*z*pow(a,3) - 1.*pow(a,4) + I*a*(5. - 4.*r)*z*pow(r,2) + 3.*(-1. + r)*pow(r,3))*pow(r + I*a*z,-2);
	}else if(nt == 1&& nr == 0&& nz == 0&& np == 1){
		hab = -0.5*a*(pow(a,2) + pow(r,2))*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
	}else if(nt == 1&& nr == 0&& nz == 1&& np == 0){
		hab = 0.;
	}else if(nt == 1&& nr == 1&& nz == 0&& np == 0){
		hab = 0.5*((-2. + r)*r + pow(a,2))*(pow(a,2) + pow(r,2))*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
	}else if(nt == 2&& nr == 0&& nz == 0&& np == 0){
		hab = -0.25*pow(r - I*a*z,2)*pow(r + I*a*z,-2)*pow(pow(a,2) + pow(r,2),2);
	}

	return hab;
}



// metric components decomposed onto tetrad basis
// Complex metric_coefficient_ORG_11(int nt, int nr, int nz, int np, double a, double r, double z){
// 	Complex hab = 0.;
// 	if(nt == 0 && nr == 0 && nz == 0 && np == 0){
// 		hab = (r - I*a*z)*(r + I*a*z);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 1){
// 		hab = (-1.*r + I*a*z)*(a + I*r*z)*pow(-1. + pow(z,2),-1);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 2){
// 		hab = -0.5*pow(r - I*a*z,2)*pow(-1. + pow(z,2),-1);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 0){
// 		hab = (r - I*a*z)*(2.*r*z - I*a*(1. + pow(z,2)));
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 1){
// 		hab = (-1.*r + I*a*z)*(I*r + a*z);
// 	}else if(nt == 0 && nr == 0 && nz == 2 && np == 0){
// 		hab = 0.5*(-1. + pow(z,2))*pow(r - I*a*z,2);
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 1){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 2 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 0){
// 		hab = a*(r - I*a*z)*(a + 2.*I*r*z + a*pow(z,2));
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 1){
// 		hab = a*pow(r - I*a*z,2);
// 	}else if(nt == 1 && nr == 0 && nz == 1 && np == 0){
// 		hab = a*(r - I*a*z)*(I*r + a*z)*(-1. + pow(z,2));
// 	}else if(nt == 1 && nr == 1 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 2 && nr == 0 && nz == 0 && np == 0){
// 		hab = -0.5*pow(a,2)*(-1. + pow(z,2))*pow(r - I*a*z,2);
// 	}
// 	return hab;
// }

// Complex metric_coefficient_ORG_13(int nt, int nr, int nz, int np, double a, double r, double z){
// 	Complex hab = 0.;
// 	if(nt == 0 && nr == 0 && nz == 0 && np == 0){
// 		hab = 2.*z*(r - I*a*z)*((-2. + r)*pow(a,2) + pow(r,3))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 1){
// 		hab = (-1.*r + I*a*z)*(I*r*(-1.*a + r)*(a + r) + a*z*(pow(a,2) + pow(r,2)) + 2.*I*(-1. + r)*pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 2){
// 		hab = 0.5*I*a*(r - I*a*z)*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 0){
// 		hab = (r - I*a*z)*(-1. + pow(z,2))*(pow(r,3) + pow(a,2)*(-1.*r + 2.*(-1. + r)*pow(z,2)))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 1){
// 		hab = 0.5*a*(-1.*r + I*a*z)*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 2 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 0){
// 		hab = -1.*z*(-1.*r + I*a*z)*((-2. + r)*r + pow(a,2))*(pow(a,2) + pow(r,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 1){
// 		hab = 0.5*I*(-1.*r + I*a*z)*((-2. + r)*r + pow(a,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 1 && nz == 1 && np == 0){
// 		hab = 0.5*(r - I*a*z)*((-2. + r)*r + pow(a,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 2 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 0){
// 		hab = -1.*(r - I*a*z)*pow(r + I*a*z,-2)*((-I)*a*(-1. + pow(z,2))*(pow(r,3) + pow(a,2)*(-1.*r + 2.*(-1. + r)*pow(z,2))) + z*pow(pow(a,2) + pow(r,2),2))*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 1){
// 		hab = 0.5*I*(r - I*a*z)*(2.*pow(a,2)*pow(r,2) + pow(r,4) - 1.*pow(a,4)*(-2. + pow(z,2))*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 1 && nr == 0 && nz == 1 && np == 0){
// 		hab = 0.5*(-1.*r + I*a*z)*(pow(a,2) + pow(r,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 1 && nr == 1 && nz == 0 && np == 0){
// 		hab = 0.5*I*a*(r - I*a*z)*((-2. + r)*r + pow(a,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 2 && nr == 0 && nz == 0 && np == 0){
// 		hab = Complex(0.,-0.5)*a*(r - I*a*z)*(pow(a,2) + pow(r,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}
// 	return hab;
// }

// Complex metric_coefficient_ORG_33(int nt, int nr, int nz, int np, double a, double r, double z){
// 	Complex hab = 0.;
// 	if(nt == 0 && nr == 0 && nz == 0 && np == 0){
// 		hab = (r - I*a*z)*(2.*r - 2.*pow(a,2) + r*pow(a,2) + I*a*z*(2. + 3.*(-2. + r)*r + pow(a,2)) - 1.*pow(r,3))*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 1){
// 		hab = 0.5*a*(-1.*r + I*a*z)*(r + Complex(0.,3.)*a*(-1. + r)*z + pow(a,2) - 2.*pow(r,2))*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 2){
// 		hab = -0.25*pow(a,2)*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 1){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 0 && nz == 2 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 0){
// 		hab = 0.5*(r - I*a*z)*((-2. + r)*r + pow(a,2))*((2. - 3.*r)*r + 4.*I*a*(-1. + r)*z + pow(a,2))*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 1){
// 		hab = 0.5*(a*(-2. + r)*r + pow(a,3))*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 1 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 2 && nz == 0 && np == 0){
// 		hab = -0.25*pow(r - I*a*z,2)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),2);
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 0){
// 		hab = 0.5*(r - I*a*z)*(r*(-1. + 2.*r)*pow(a,2) + I*(3. - 4.*r)*z*pow(a,3) - 1.*pow(a,4) + I*a*(5. - 4.*r)*z*pow(r,2) + 3.*(-1. + r)*pow(r,3))*pow(r + I*a*z,-2);
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 1){
// 		hab = -0.5*a*(pow(a,2) + pow(r,2))*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
// 	}else if(nt == 1 && nr == 0 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 1 && nr == 1 && nz == 0 && np == 0){
// 		hab = 0.5*((-2. + r)*r + pow(a,2))*(pow(a,2) + pow(r,2))*pow(r - I*a*z,2)*pow(r + I*a*z,-2);
// 	}else if(nt == 2 && nr == 0 && nz == 0 && np == 0){
// 		hab = -0.25*pow(r - I*a*z,2)*pow(r + I*a*z,-2)*pow(pow(a,2) + pow(r,2),2);
// 	}
// 	return hab;
// }

// Complex metric_coefficient_IRG_22(int nt, int nr, int nz, int np, double a, double r, double z){
// 	Complex hab = 0.;
// 	if(nt == 0 && nr == 0 && nz == 0 && np == 0){
// 		hab = pow(pow(r,2) + pow(a,2)*pow(z,2),-1);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 1){
// 		hab = (a + I*r*z)*pow(r - I*a*z,-1)*pow(r + I*a*z,-2)*pow(-1. + pow(z,2),-1);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 2){
// 		hab = -0.5*pow(r + I*a*z,-2)*pow(-1. + pow(z,2),-1);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 0){
// 		hab = (2.*r*z - I*a*(1. + pow(z,2)))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 1){
// 		hab = (I*r + a*z)*pow(r - I*a*z,-1)*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 0 && nz == 2 && np == 0){
// 		hab = 0.5*(-1. + pow(z,2))*pow(r + I*a*z,-2);
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 1){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 2 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 0){
// 		hab = -1.*a*(a + 2.*I*r*z + a*pow(z,2))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2);
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 1){
// 		hab = a*pow(r + I*a*z,-2);
// 	}else if(nt == 1 && nr == 0 && nz == 1 && np == 0){
// 		hab = -1.*a*(I*r + a*z)*(-1. + pow(z,2))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2);
// 	}else if(nt == 1 && nr == 1 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 2 && nr == 0 && nz == 0 && np == 0){
// 		hab = -0.5*pow(a,2)*(-1. + pow(z,2))*pow(r + I*a*z,-2);
// 	}
// 	return hab;
// }

// Complex metric_coefficient_IRG_24(int nt, int nr, int nz, int np, double a, double r, double z){
// 	Complex hab = 0.;
// 	if(nt == 0 && nr == 0 && nz == 0 && np == 0){
// 		hab = -4.*r*z*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 1){
// 		hab = 2.*((-I)*r*((-2. + r)*r + pow(a,2)) + a*z*(pow(a,2) + pow(r,2)))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 2){
// 		hab = I*a*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 0){
// 		hab = -1.*r*pow(r - I*a*z,-1)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 1){
// 		hab = a*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 0 && nz == 2 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 0){
// 		hab = 1.4142135623730951*z*(pow(a,2) + pow(r,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow(1. - 1.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 1){
// 		hab = I*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 0 && nr == 1 && nz == 1 && np == 0){
// 		hab = 0.5*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),0.5);
// 	}else if(nt == 0 && nr == 2 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 0){
// 		hab = 1.4142135623730951*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*(I*a*r*((-2. + r)*r + pow(a,2))*(-1. + pow(z,2)) + z*pow(pow(a,2) + pow(r,2),2))*pow(1. - 1.*pow(z,2),-0.5);
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 1){
// 		hab = I*(2.*pow(a,2)*pow(r,2) + pow(r,4) - 1.*pow(a,4)*(-2. + pow(z,2))*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 1 && nr == 0 && nz == 1 && np == 0){
// 		hab = (pow(a,2) + pow(r,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),-0.5);
// 	}else if(nt == 1 && nr == 1 && nz == 0 && np == 0){
// 		hab = Complex(0.,-0.5)*a*(pow(r,2) + pow(a,2)*pow(z,2))*pow(r - I*a*z,-1)*pow(r + I*a*z,-2)*pow(2. - 2.*pow(z,2),0.5);
// 	}else if(nt == 2 && nr == 0 && nz == 0 && np == 0){
// 		hab = (-I)*a*(pow(a,2) + pow(r,2))*(-1. + pow(z,2))*(pow(r,2) + pow(a,2)*pow(z,2))*pow(-1.*r + I*a*z,-1)*pow(r + I*a*z,-2)*pow((-2. + r)*r + pow(a,2),-1)*pow(2. - 2.*pow(z,2),-0.5);
// 	}
// 	return hab;
// }

// Complex metric_coefficient_IRG_44(int nt, int nr, int nz, int np, double a, double r, double z){
// 	Complex hab = 0.;
// 	if(nt == 0 && nr == 0 && nz == 0 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 1){
// 		hab = -2.*a*(r*(-3. + 2.*r) - I*a*(-1. + r)*z + pow(a,2))*pow(-1.*r + I*a*z,-1)*pow((-2. + r)*r + pow(a,2),-2);
// 	}else if(nt == 0 && nr == 0 && nz == 0 && np == 2){
// 		hab = -1.*pow(a,2)*pow((-2. + r)*r + pow(a,2),-2);
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 0 && nz == 1 && np == 1){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 0 && nz == 2 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 0){
// 		hab = 2.*pow(r - I*a*z,-1);
// 	}else if(nt == 0 && nr == 1 && nz == 0 && np == 1){
// 		hab = -2.*a*pow((-2. + r)*r + pow(a,2),-1);
// 	}else if(nt == 0 && nr == 1 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 0 && nr == 2 && nz == 0 && np == 0){
// 		hab = -1.;
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 0){
// 		hab = 2.*(-1.*pow(a,2) + pow(r,2) - 1.*((-2. + r)*r + pow(a,2))*(pow(a,2) + pow(r,2))*pow(-1.*r + I*a*z,-1))*pow((-2. + r)*r + pow(a,2),-2);
// 	}else if(nt == 1 && nr == 0 && nz == 0 && np == 1){
// 		hab = -2.*a*(pow(a,2) + pow(r,2))*pow((-2. + r)*r + pow(a,2),-2);
// 	}else if(nt == 1 && nr == 0 && nz == 1 && np == 0){
// 		hab = 0.;
// 	}else if(nt == 1 && nr == 1 && nz == 0 && np == 0){
// 		hab = -2. - 4.*r*pow((-2. + r)*r + pow(a,2),-1);
// 	}else if(nt == 2 && nr == 0 && nz == 0 && np == 0){
// 		hab = -1.*pow((-2. + r)*r + pow(a,2),-2)*pow(pow(a,2) + pow(r,2),2);
// 	}
// 	return hab;
// }



// // metric coefficients when decomposed onto tetrad
// ComplexTensor metric_coefficients_11(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*(r[jr] - I*a*z[jz])*(I*r[jr]*z[jz] + a*(-1. + 2.*pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*pow(r[jr] - I*a*z[jz],2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-1.*r[jr] + I*a*z[jz])*(r[jr]*z[jz] + I*(a - 2.*a*pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-1.*r[jr] + I*a*z[jz])*(I*r[jr] + a*z[jz])*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*pow(r[jr] - I*a*z[jz],2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (r[jr] - I*a*z[jz])*pow(a,2);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*a*pow(r[jr] - I*a*z[jz],2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(r[jr] - I*a*z[jz])*(I*r[jr] + a*z[jz])*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*pow(a,2)*pow(r[jr] - I*a*z[jz],2);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_11_dz(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (2.*a*r[jr]*z[jz]*(1. + pow(z[jz],2)) + I*pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) + 4.*a*r[jr]*pow(z[jz],3) + I*pow(a,2)*(1. - 3.*pow(z[jz],2) - 2.*pow(z[jz],4)))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*(r[jr] - I*a*z[jz])*(2.*r[jr]*z[jz] - I*a*(1. + pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (Complex(0.,-2.)*a*r[jr]*z[jz]*(1. + 3.*pow(z[jz],2)) + pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) - 1.*pow(a,2)*(-1. + 3.*pow(z[jz],2) + 2.*pow(z[jz],4)))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (4.*I*z[jz]*pow(r[jr],2) - 2.*I*z[jz]*pow(a,2)*(1. + pow(z[jz],2)) + 2.*a*r[jr]*(1. + 3.*pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (r[jr] - I*a*z[jz])*(2.*r[jr]*z[jz] - I*a*(1. + pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-I)*pow(a,3);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*a*(r[jr] - I*a*z[jz])*((-I)*a + r[jr]*z[jz])*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(2.*I*z[jz]*pow(a,2) - 2.*I*z[jz]*pow(r[jr],2) - 2.*a*r[jr]*(1. + pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*(-1.*r[jr] + I*a*z[jz])*pow(a,3);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_11_dz2(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(-6.*a*r[jr]*pow(z[jz],2)*(1. + pow(z[jz],2)) + 2.*I*pow(a,2)*(5. + pow(z[jz],2))*pow(z[jz],3) - 6.*I*pow(r[jr],2)*(z[jz] + pow(z[jz],3)) - 1.*a*r[jr]*(1. + 8.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (Complex(0.,-12.)*a*r[jr]*z[jz]*(1. + pow(z[jz],2)) + 2.*pow(r[jr],2)*(1. + 5.*pow(z[jz],2)) - 1.*pow(a,2)*(1. + 8.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(-6.*z[jz]*pow(r[jr],2)*(1. + pow(z[jz],2)) + 2.*pow(a,2)*(5. + pow(z[jz],2))*pow(z[jz],3) + I*a*r[jr]*(1. + 14.*pow(z[jz],2) + 9.*pow(z[jz],4)))*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(-12.*a*r[jr]*z[jz]*(1. + pow(z[jz],2)) - 2.*I*pow(r[jr],2)*(1. + 5.*pow(z[jz],2)) + I*pow(a,2)*(1. + 8.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (Complex(0.,12.)*a*r[jr]*z[jz]*(1. + pow(z[jz],2)) - 2.*pow(r[jr],2)*(1. + 5.*pow(z[jz],2)) + pow(a,2)*(1. + 8.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*a*(Complex(0.,-2.)*a*r[jr]*z[jz]*(3. + pow(z[jz],2)) - 1.*pow(a,2)*(1. + 3.*pow(z[jz],2)) + pow(r[jr],2)*(1. + 3.*pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*a*(-2.*a*r[jr]*z[jz]*(3. + pow(z[jz],2)) + I*(a - 1.*r[jr])*(a + r[jr])*(1. + 3.*pow(z[jz],2)))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*pow(a,4);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_13(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.8284271247461903*(-1. + r[jr])*z[jz]*(r[jr] - I*a*z[jz])*pow(a,2)*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-0.5);

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (r[jr] - I*a*z[jz])*(I*r[jr]*(-1.*a + r[jr])*(a + r[jr]) - 1.*z[jz]*pow(a,3)*(-1. + pow(z[jz],2)) + 2.*I*(-1. + r[jr])*pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,-0.35355339059327373)*a*(-1.*r[jr] + I*a*z[jz])*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (r[jr] - I*a*z[jz])*(pow(r[jr],3) + pow(a,2)*(-1.*r[jr] + 2.*(-1. + r[jr])*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.35355339059327373*a*(-1.*r[jr] + I*a*z[jz])*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = z[jz]*(r[jr] - I*a*z[jz])*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,0.35355339059327373)*(-1.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.35355339059327373*(-1.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(-1.*r[jr] + I*a*z[jz])*(I*r[jr]*(-1.*a + r[jr])*(a + r[jr]) + a*z[jz]*(pow(a,2) + pow(r[jr],2)) + 2.*I*(-1. + r[jr])*pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,-0.5)*(r[jr] - I*a*z[jz])*(2.*pow(a,2)*pow(r[jr],2) + pow(r[jr],4) - 1.*pow(a,4)*(-2. + pow(z[jz],2))*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.35355339059327373*(-1.*r[jr] + I*a*z[jz])*(pow(a,2) + pow(r[jr],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*a*(-1.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,-0.5)*a*(-1.*r[jr] + I*a*z[jz])*(pow(a,2) + pow(r[jr],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_13_dz(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.8284271247461903*(-1. + r[jr])*pow(a,2)*(-1.*pow(r[jr],2) - Complex(0.,3.)*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (Complex(0.,-3.)*z[jz]*pow(r[jr],5) + 3.*a*pow(r[jr],4)*(-1. + pow(z[jz],2)) + I*z[jz]*pow(a,2)*pow(r[jr],2)*(4. + 2.*(1. - 3.*r[jr])*pow(z[jz],2)) - 1.*r[jr]*pow(a,3)*(-1. + pow(z[jz],2))*(2.*r[jr] - 6.*(-1. + r[jr])*pow(z[jz],2)) + pow(a,5)*(-1. + pow(z[jz],2))*pow(z[jz],4) + I*z[jz]*pow(a,4)*(2.*r[jr] - 2.*(-1. + 2.*r[jr])*pow(z[jz],2) + (-4.*(-1. + r[jr]) + 3.*r[jr])*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*a*(Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) + z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2) - 1.*(pow(a,2) + 4.*pow(r[jr],2))*pow(z[jz],2) - 2.*pow(a,2)*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (Complex(0.,-3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],3) + pow(a,2)*(-1.*r[jr] + 2.*(-1. + r[jr])*pow(z[jz],2))) + z[jz]*pow(r[jr],2)*(-3.*pow(r[jr],3) - 1.*pow(a,2)*(-4. + r[jr] + 2.*(-1. + r[jr])*pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(r[jr],3)*(-1. + 4.*pow(z[jz],2)) + pow(a,2)*(r[jr] - 2.*(1. + r[jr])*pow(z[jz],2) + 4.*(-1. + r[jr])*pow(z[jz],4))))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*a*(Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) + z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2) - 1.*(pow(a,2) + 4.*pow(r[jr],2))*pow(z[jz],2) - 2.*pow(a,2)*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(-1.*pow(r[jr],2) - Complex(0.,3.)*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*((-2. + r[jr])*r[jr] + pow(a,2))*(-1.*z[jz]*pow(a,2)*pow(r[jr],2)*(-1. + 4.*pow(z[jz],2)) - Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*pow(a,4)*(pow(z[jz],3) + 2.*pow(z[jz],5)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*((-2. + r[jr])*r[jr] + pow(a,2))*(-1.*z[jz]*pow(a,2)*pow(r[jr],2)*(-1. + 4.*pow(z[jz],2)) - Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*pow(a,4)*(pow(z[jz],3) + 2.*pow(z[jz],5)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(I*z[jz]*pow(r[jr],5) + a*pow(r[jr],4)*(1. - 3.*(-1. + pow(z[jz],2))) + I*z[jz]*pow(a,4)*(-2.*r[jr] + (-2. + 3.*r[jr])*pow(z[jz],2)) + I*z[jz]*pow(a,2)*pow(r[jr],2)*(-4. - 1.*r[jr] + (2. + 3.*r[jr])*pow(z[jz],2)) + pow(a,5)*pow(z[jz],4) + r[jr]*pow(a,3)*(r[jr] + 6.*(-1. + pow(z[jz],2))*pow(z[jz],2) - 1.*r[jr]*(3. - 9.*pow(z[jz],2) + 5.*pow(z[jz],4))))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*(Complex(0.,-3.)*a*r[jr]*(-1.*pow(r[jr],2) + pow(a,2)*(-2. + pow(z[jz],2)))*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(2.*pow(a,2)*pow(r[jr],2)*(1. - 4.*pow(z[jz],2)) + pow(r[jr],4)*(1. - 4.*pow(z[jz],2)) - 1.*pow(a,4)*pow(z[jz],2)*(2. + pow(z[jz],2))) + z[jz]*pow(r[jr],2)*(6.*pow(a,2)*pow(r[jr],2) + 3.*pow(r[jr],4) + pow(a,4)*(4. - 2.*pow(z[jz],2) + pow(z[jz],4))))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*(pow(a,2) + pow(r[jr],2))*(-1.*z[jz]*pow(a,2)*pow(r[jr],2)*(-1. + 4.*pow(z[jz],2)) - Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*pow(a,4)*(pow(z[jz],3) + 2.*pow(z[jz],5)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*a*((-2. + r[jr])*r[jr] + pow(a,2))*(z[jz]*pow(r[jr],2)*(pow(r[jr],2) - 1.*pow(a,2)*(-2. + pow(z[jz],2))) + Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2) - 1.*(pow(a,2) + 2.*pow(r[jr],2))*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*a*(pow(a,2) + pow(r[jr],2))*(-1.*z[jz]*pow(r[jr],4) + z[jz]*pow(a,2)*pow(r[jr],2)*(-2. + pow(z[jz],2)) - Complex(0.,3.)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(a,2)*pow(z[jz],2) + pow(r[jr],2)*(-1. + 2.*pow(z[jz],2))))*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_13_dz2(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.8284271247461903*(-1. + r[jr])*pow(a,2)*(-3.*z[jz]*pow(r[jr],3) + Complex(0.,3.)*a*pow(r[jr],2)*(2. - 3.*pow(z[jz],2)) - 1.*r[jr]*z[jz]*pow(a,2)*(-6. + 19.*pow(z[jz],2) - 10.*pow(z[jz],4)) - I*pow(a,3)*(1. + 2.*pow(z[jz],2))*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-2.5);

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-3.*z[jz]*pow(a,3)*pow(r[jr],3)*(-1. + pow(z[jz],2)) + I*pow(r[jr],2)*(3.*pow(r[jr],4)*(1. + 4.*pow(z[jz],2)) - 3.*pow(a,4)*(2. - 5.*pow(z[jz],2) + 3.*pow(z[jz],4)) + r[jr]*pow(a,2)*(-4. + r[jr] + 2.*(-11. + 5.*r[jr])*pow(z[jz],2) + 4.*(-1. + r[jr])*pow(z[jz],4))) - 1.*a*r[jr]*z[jz]*(15.*pow(r[jr],4)*(-1. + 2.*pow(z[jz],2)) + r[jr]*pow(a,2)*(-5.*(-4. + r[jr]) + 2.*(-17. + 2.*r[jr])*pow(z[jz],2) + 16.*(-1. + r[jr])*pow(z[jz],4)) + pow(a,4)*(6. - 25.*pow(z[jz],2) + 29.*pow(z[jz],4) - 10.*pow(z[jz],6))) + z[jz]*pow(a,3)*(pow(r[jr],3)*(-2. + 7.*pow(z[jz],2) - 20.*pow(z[jz],4)) + pow(a,2)*(2.*r[jr] - 7.*r[jr]*pow(z[jz],2) + 2.*(9. + r[jr])*pow(z[jz],4) - 12.*(-1. + r[jr])*pow(z[jz],6))) - I*pow(a,2)*(5.*pow(r[jr],4)*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4)) + pow(a,4)*pow(z[jz],4)*(-1. - 1.*pow(z[jz],2) + 2.*pow(z[jz],4)) + r[jr]*pow(a,2)*(-10.*r[jr] + 35.*r[jr]*pow(z[jz],2) + 2.*(27. - 32.*r[jr])*pow(z[jz],4) + 24.*(-1. + r[jr])*pow(z[jz],6))))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*a*(Complex(0.,15.)*a*z[jz]*pow(r[jr],4)*(1. - 2.*pow(z[jz],2)) - 3.*pow(r[jr],5)*(1. + 4.*pow(z[jz],2)) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 10.*pow(z[jz],2) - 28.*pow(z[jz],4)) - 3.*r[jr]*pow(a,4)*(9. - 4.*pow(z[jz],2))*pow(z[jz],4) - 1.*pow(a,2)*pow(r[jr],3)*(2. + 11.*pow(z[jz],2) + 2.*pow(z[jz],4) - 5.*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4))) - Complex(0.,3.)*pow(a,5)*(3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (3.*pow(r[jr],6)*(1. + 4.*pow(z[jz],2)) + pow(a,2)*pow(r[jr],3)*(-4. + r[jr] + 2.*(-11. + 5.*r[jr])*pow(z[jz],2) + 4.*(-1. + r[jr])*pow(z[jz],4)) + I*a*z[jz]*pow(r[jr],2)*(15.*pow(r[jr],3)*(-1. + 2.*pow(z[jz],2)) + pow(a,2)*(-5.*(-4. + r[jr]) + 2.*(-17. + 2.*r[jr])*pow(z[jz],2) + 16.*(-1. + r[jr])*pow(z[jz],4))) - I*z[jz]*pow(a,3)*(pow(r[jr],3)*(-2. + 7.*pow(z[jz],2) - 20.*pow(z[jz],4)) + pow(a,2)*(2.*r[jr] - 7.*r[jr]*pow(z[jz],2) + 2.*(9. + r[jr])*pow(z[jz],4) - 12.*(-1. + r[jr])*pow(z[jz],6))) - 1.*r[jr]*pow(a,2)*(5.*pow(r[jr],3)*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4)) + pow(a,2)*(-10.*r[jr] + 35.*r[jr]*pow(z[jz],2) + 2.*(27. - 32.*r[jr])*pow(z[jz],4) + 24.*(-1. + r[jr])*pow(z[jz],6))))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*a*(Complex(0.,15.)*a*z[jz]*pow(r[jr],4)*(1. - 2.*pow(z[jz],2)) - 3.*pow(r[jr],5)*(1. + 4.*pow(z[jz],2)) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 10.*pow(z[jz],2) - 28.*pow(z[jz],4)) - 3.*r[jr]*pow(a,4)*(9. - 4.*pow(z[jz],2))*pow(z[jz],4) - 1.*pow(a,2)*pow(r[jr],3)*(2. + 11.*pow(z[jz],2) + 2.*pow(z[jz],4) - 5.*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4))) - Complex(0.,3.)*pow(a,5)*(3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(-3.*z[jz]*pow(r[jr],3) + Complex(0.,3.)*a*pow(r[jr],2)*(2. - 3.*pow(z[jz],2)) - 1.*r[jr]*z[jz]*pow(a,2)*(-6. + 19.*pow(z[jz],2) - 10.*pow(z[jz],4)) - I*pow(a,3)*(1. + 2.*pow(z[jz],2))*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,-0.5)*((-2. + r[jr])*r[jr] + pow(a,2))*(Complex(0.,15.)*a*z[jz]*pow(r[jr],4)*(1. - 2.*pow(z[jz],2)) - 3.*pow(r[jr],5)*(1. + 4.*pow(z[jz],2)) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 10.*pow(z[jz],2) - 28.*pow(z[jz],4)) - 3.*r[jr]*pow(a,4)*(9. - 4.*pow(z[jz],2))*pow(z[jz],4) - 1.*pow(a,2)*pow(r[jr],3)*(2. + 11.*pow(z[jz],2) + 2.*pow(z[jz],4) - 5.*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4))) - Complex(0.,3.)*pow(a,5)*(3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*((-2. + r[jr])*r[jr] + pow(a,2))*(Complex(0.,15.)*a*z[jz]*pow(r[jr],4)*(1. - 2.*pow(z[jz],2)) - 3.*pow(r[jr],5)*(1. + 4.*pow(z[jz],2)) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 10.*pow(z[jz],2) - 28.*pow(z[jz],4)) - 3.*r[jr]*pow(a,4)*(9. - 4.*pow(z[jz],2))*pow(z[jz],4) - 1.*pow(a,2)*pow(r[jr],3)*(2. + 11.*pow(z[jz],2) + 2.*pow(z[jz],4) - 5.*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4))) - Complex(0.,3.)*pow(a,5)*(3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(-3.*a*z[jz]*(pow(a,2) + pow(r[jr],2))*pow(r[jr],3) - I*pow(r[jr],2)*(pow(r[jr],4)*(1. + 2.*pow(z[jz],2)) + pow(a,4)*(-6. + 9.*pow(z[jz],2)) + r[jr]*pow(a,2)*(-4. - 3.*r[jr] + (-2. + 9.*r[jr])*pow(z[jz],2))) + z[jz]*pow(a,3)*(pow(a,2)*(-2.*r[jr] + 5.*r[jr]*pow(z[jz],2) - 6.*pow(z[jz],4)) + pow(r[jr],3)*(2. - 5.*pow(z[jz],2) + 6.*pow(z[jz],4))) - 1.*a*r[jr]*z[jz]*(pow(r[jr],4)*(5. - 8.*pow(z[jz],2)) + pow(a,4)*(-6. + 19.*pow(z[jz],2) - 10.*pow(z[jz],4)) + r[jr]*pow(a,2)*(-20. + 9.*r[jr] + (38. - 11.*r[jr])*pow(z[jz],2) + 2.*(-6. + r[jr])*pow(z[jz],4))) - I*pow(a,2)*(pow(r[jr],4)*(-10. + 25.*pow(z[jz],2) - 12.*pow(z[jz],4)) + pow(a,4)*(pow(z[jz],4) + 2.*pow(z[jz],6)) + r[jr]*pow(a,2)*(10.*r[jr] - 25.*r[jr]*pow(z[jz],2) + (-18. + 31.*r[jr])*pow(z[jz],4) + 2.*(6. - 5.*r[jr])*pow(z[jz],6))))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*(Complex(0.,15.)*a*z[jz]*pow(r[jr],6)*(1. - 2.*pow(z[jz],2)) - 3.*pow(r[jr],7)*(1. + 4.*pow(z[jz],2)) + I*z[jz]*pow(a,3)*pow(r[jr],4)*(28. - 53.*pow(z[jz],2) - 20.*pow(z[jz],4)) - 3.*r[jr]*pow(a,6)*pow(z[jz],4)*(12. - 9.*pow(z[jz],2) + 2.*pow(z[jz],4)) - 1.*pow(a,4)*pow(r[jr],3)*(4. + 10.*pow(z[jz],2) + pow(z[jz],4) - 10.*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4))) - 1.*pow(a,2)*pow(r[jr],5)*(6. + 24.*pow(z[jz],2) - 5.*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4))) - Complex(0.,3.)*pow(a,7)*(4. + pow(z[jz],2))*pow(z[jz],5) + I*z[jz]*pow(a,5)*pow(r[jr],2)*(20. - 46.*pow(z[jz],2) + 17.*pow(z[jz],4) - 2.*(2. - 7.*pow(z[jz],2) + 20.*pow(z[jz],4)) - 6.*pow(z[jz],6)))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*(pow(a,2) + pow(r[jr],2))*(Complex(0.,15.)*a*z[jz]*pow(r[jr],4)*(1. - 2.*pow(z[jz],2)) - 3.*pow(r[jr],5)*(1. + 4.*pow(z[jz],2)) + I*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 10.*pow(z[jz],2) - 28.*pow(z[jz],4)) - 3.*r[jr]*pow(a,4)*(9. - 4.*pow(z[jz],2))*pow(z[jz],4) - 1.*pow(a,2)*pow(r[jr],3)*(2. + 11.*pow(z[jz],2) + 2.*pow(z[jz],4) - 5.*(2. - 7.*pow(z[jz],2) + 2.*pow(z[jz],4))) - Complex(0.,3.)*pow(a,5)*(3. + 2.*pow(z[jz],2))*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,-0.5)*a*((-2. + r[jr])*r[jr] + pow(a,2))*((-I)*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 14.*pow(z[jz],2)) + pow(r[jr],5)*(1. + 2.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],4)*(-5. + 8.*pow(z[jz],2)) + pow(a,2)*pow(r[jr],3)*(-8. + 26.*pow(z[jz],2) - 12.*pow(z[jz],4)) - 3.*r[jr]*pow(a,4)*(-3. + 2.*pow(z[jz],2))*pow(z[jz],4) + Complex(0.,3.)*pow(a,5)*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*I*a*(pow(a,2) + pow(r[jr],2))*((-I)*z[jz]*pow(a,3)*pow(r[jr],2)*(8. - 14.*pow(z[jz],2)) + pow(r[jr],5)*(1. + 2.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],4)*(-5. + 8.*pow(z[jz],2)) + pow(a,2)*pow(r[jr],3)*(-8. + 26.*pow(z[jz],2) - 12.*pow(z[jz],4)) - 3.*r[jr]*pow(a,4)*(-3. + 2.*pow(z[jz],2))*pow(z[jz],4) + Complex(0.,3.)*pow(a,5)*pow(z[jz],5))*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_33(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (r[jr] - I*a*z[jz])*(-2.*r[jr] + 2.*pow(a,2) - 1.*r[jr]*pow(a,2) - I*a*z[jz]*(2. + 3.*(-2. + r[jr])*r[jr] + pow(a,2)) + pow(r[jr],3))*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*a*(r[jr] - I*a*z[jz])*(r[jr] + Complex(0.,3.)*a*(-1. + r[jr])*z[jz] + pow(a,2) - 2.*pow(r[jr],2))*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.25*pow(a,2)*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*(r[jr] - I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*((2. - 3.*r[jr])*r[jr] + 4.*I*a*(-1. + r[jr])*z[jz] + pow(a,2))*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*(a*(-2. + r[jr])*r[jr] + pow(a,3))*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.25*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*(r[jr] - I*a*z[jz])*(r[jr]*(-1. + 2.*r[jr])*pow(a,2) + I*(3. - 4.*r[jr])*z[jz]*pow(a,3) - 1.*pow(a,4) + I*a*(5. - 4.*r[jr])*z[jz]*pow(r[jr],2) + 3.*(-1. + r[jr])*pow(r[jr],3))*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*a*(pow(a,2) + pow(r[jr],2))*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2))*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.25*pow(r[jr] - I*a*z[jz],2)*pow(r[jr] + I*a*z[jz],-2)*pow(pow(a,2) + pow(r[jr],2),2)*pow(-1. + pow(z[jz],2),-1);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_33_dz(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*(z[jz]*pow(r[jr],2)*(2.*pow(a,2) - 1.*r[jr]*(2. + pow(a,2)) + pow(r[jr],3)) + I*a*r[jr]*(pow(a,2)*(-3. + r[jr] + 3.*pow(z[jz],2) - 2.*r[jr]*pow(z[jz],2)) + r[jr]*(2. - 4.*pow(z[jz],2) + 3.*r[jr]*(1. - 1.*r[jr] + pow(z[jz],2)))) - 1.*z[jz]*pow(a,2)*(pow(a,2)*(1. + r[jr] - 3.*pow(z[jz],2)) + r[jr]*(2. + r[jr]*(-9. + 5.*r[jr] + 9.*pow(z[jz],2) - 6.*r[jr]*pow(z[jz],2)))) - I*(2. + 3.*(-2. + r[jr])*r[jr] + pow(a,2))*pow(a,3)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*a*(2.*z[jz]*(r[jr]*(-1. + 2.*r[jr]) - 1.*pow(a,2))*pow(r[jr],2) + Complex(0.,3.)*a*r[jr]*(r[jr]*(2. + r[jr]*(-3. + pow(z[jz],2))) - 1.*pow(a,2)*(-1. + pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(a,2)*(-1. + 3.*pow(z[jz],2)) + r[jr]*(-10. + 11.*r[jr] + 3.*(4. - 5.*r[jr])*pow(z[jz],2))) - 6.*I*(-1. + r[jr])*pow(a,3)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*(r[jr] - I*a*z[jz])*pow(a,2)*(z[jz]*pow(r[jr],2) + 2.*I*a*r[jr]*(-1. + pow(z[jz],2)) + pow(a,2)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*((-2. + r[jr])*r[jr] + pow(a,2))*(2.*z[jz]*(r[jr]*(-2. + 3.*r[jr]) - 1.*pow(a,2))*pow(r[jr],2) - 1.*z[jz]*pow(a,2)*(pow(a,2)*(-1. + 3.*pow(z[jz],2)) + r[jr]*(-14. + 15.*r[jr] + 3.*(6. - 7.*r[jr])*pow(z[jz],2))) + I*a*r[jr]*(-3.*pow(a,2)*(-1. + pow(z[jz],2)) + r[jr]*(10. - 13.*r[jr] + (-2. + 5.*r[jr])*pow(z[jz],2))) - 8.*I*(-1. + r[jr])*pow(a,3)*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(-1.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(-1.*z[jz]*pow(r[jr],2) - 2.*I*a*r[jr]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*(-1.*r[jr] + I*a*z[jz])*(-1.*z[jz]*pow(r[jr],2) - 2.*I*a*r[jr]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*(2.*z[jz]*pow(r[jr],2)*(r[jr]*(-1. + 2.*r[jr])*pow(a,2) - 1.*pow(a,4) + 3.*(-1. + r[jr])*pow(r[jr],3)) + I*a*r[jr]*(2.*r[jr]*pow(a,2)*(3. + r[jr]*(-5. + pow(z[jz],2))) - 3.*pow(a,4)*(-1. + pow(z[jz],2)) + pow(r[jr],3)*(14. - 13.*r[jr] + (-4. + 5.*r[jr])*pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(a,4)*(-1. + 3.*pow(z[jz],2)) - 2.*r[jr]*pow(a,2)*(5. - 7.*r[jr] + 3.*(-2. + 3.*r[jr])*pow(z[jz],2)) - 3.*pow(r[jr],3)*(6. - 5.*r[jr] + (-8. + 7.*r[jr])*pow(z[jz],2))) - 2.*I*pow(a,3)*((-3. + 4.*r[jr])*pow(a,2) + (-5. + 4.*r[jr])*pow(r[jr],2))*pow(z[jz],4))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*a*(-1.*r[jr] + I*a*z[jz])*(pow(a,2) + pow(r[jr],2))*(-1.*z[jz]*pow(r[jr],2) - 2.*I*a*r[jr]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-1.*r[jr] + I*a*z[jz])*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2))*(-1.*z[jz]*pow(r[jr],2) - 2.*I*a*r[jr]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*(-1.*r[jr] + I*a*z[jz])*(-1.*z[jz]*pow(r[jr],2) - 2.*I*a*r[jr]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*pow(z[jz],3))*pow(r[jr] + I*a*z[jz],-3)*pow(pow(a,2) + pow(r[jr],2),2)*pow(-1. + pow(z[jz],2),-2);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_33_dz2(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-3)*((-I)*a*(3.*r[jr] - I*a*z[jz])*(r[jr] + I*a*z[jz])*(1. - 1.*pow(z[jz],2))*(2.*z[jz]*(2.*pow(a,2) - 1.*r[jr]*(2. + pow(a,2)) + pow(r[jr],3)) - I*a*(2. + 3.*(-2. + r[jr])*r[jr] + pow(a,2))*(1. + pow(z[jz],2))) - 1.*(r[jr] - I*a*z[jz])*(I*a*(2. + 3.*(-2. + r[jr])*r[jr])*z[jz]*(3. + pow(z[jz],2)) + I*z[jz]*pow(a,3)*(3. + pow(z[jz],2)) + (-2. + r[jr])*pow(a,2)*(1. + 3.*pow(z[jz],2)) - 1.*r[jr]*(-2. + pow(r[jr],2))*(1. + 3.*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],2) - 1.*(-5.*r[jr] + I*a*z[jz])*pow(a,2)*(2.*r[jr] - 2.*pow(a,2) + r[jr]*pow(a,2) + I*a*z[jz]*(2. + 3.*(-2. + r[jr])*r[jr] + pow(a,2)) - 1.*pow(r[jr],3))*pow(-1. + pow(z[jz],2),2));

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*((r[jr] + pow(a,2) - 2.*pow(r[jr],2))*pow(r[jr],3)*(1. + 3.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],2)*(pow(a,2)*(-5. + 9.*pow(z[jz],2)) + r[jr]*(-14. + 19.*r[jr] + 3.*(2. - 5.*r[jr])*pow(z[jz],2))) - 3.*(-1. + r[jr])*pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) + I*z[jz]*pow(a,3)*((10. - 11.*r[jr])*r[jr] + pow(a,2) - 3.*((14. - 15.*r[jr])*r[jr] + pow(a,2))*pow(z[jz],2) + 6.*((4. - 5.*r[jr])*r[jr] + pow(a,2))*pow(z[jz],4)) - 1.*r[jr]*pow(a,2)*(r[jr]*(14. - 19.*r[jr] + 3.*(-8. + 13.*r[jr])*pow(z[jz],2) - 6.*pow(z[jz],4)) + pow(a,2)*(5. - 15.*pow(z[jz],2) + 6.*pow(z[jz],4))))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*pow(a,2)*(8.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(r[jr],4)*(1. + 3.*pow(z[jz],2)) + pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*pow(a,2)*pow(r[jr],2)*(4. - 9.*pow(z[jz],2) + pow(z[jz],4)) + 4.*I*r[jr]*z[jz]*pow(a,3)*(1. - 4.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = ((-2. + r[jr])*r[jr] + pow(a,2))*(-1.*((2. - 3.*r[jr])*r[jr] + pow(a,2))*pow(r[jr],3)*(1. + 3.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],2)*(pow(a,2)*(5. - 9.*pow(z[jz],2)) + r[jr]*(22. - 27.*r[jr] + (-14. + 23.*r[jr])*pow(z[jz],2))) + 4.*(-1. + r[jr])*pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - I*z[jz]*pow(a,3)*(pow(a,2)*(1. - 3.*pow(z[jz],2) + 6.*pow(z[jz],4)) + r[jr]*(14. - 15.*r[jr] + (-58. + 61.*r[jr])*pow(z[jz],2) + 6.*(6. - 7.*r[jr])*pow(z[jz],4))) - 1.*r[jr]*pow(a,2)*(pow(a,2)*(-5. + 15.*pow(z[jz],2) - 6.*pow(z[jz],4)) + r[jr]*(-22. + 27.*r[jr] + 3.*(14. - 19.*r[jr])*pow(z[jz],2) + 2.*(2. + r[jr])*pow(z[jz],4))))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*(a*(-2. + r[jr])*r[jr] + pow(a,3))*(8.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(r[jr],4)*(1. + 3.*pow(z[jz],2)) + pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*pow(a,2)*pow(r[jr],2)*(4. - 9.*pow(z[jz],2) + pow(z[jz],4)) + 4.*I*r[jr]*z[jz]*pow(a,3)*(1. - 4.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*(8.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(r[jr],4)*(1. + 3.*pow(z[jz],2)) + pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*pow(a,2)*pow(r[jr],2)*(4. - 9.*pow(z[jz],2) + pow(z[jz],4)) + 4.*I*r[jr]*z[jz]*pow(a,3)*(1. - 4.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),2)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (pow(r[jr],3)*((1. - 2.*r[jr])*r[jr]*pow(a,2) + pow(a,4) - 3.*(-1. + r[jr])*pow(r[jr],3))*(1. + 3.*pow(z[jz],2)) + I*a*z[jz]*pow(r[jr],2)*(pow(a,4)*(-5. + 9.*pow(z[jz],2)) + pow(r[jr],3)*(-30. + 27.*r[jr] + 22.*pow(z[jz],2) - 23.*r[jr]*pow(z[jz],2)) + 2.*r[jr]*pow(a,2)*(-7. + 11.*r[jr] + 3.*pow(z[jz],2) - 7.*r[jr]*pow(z[jz],2))) - 1.*pow(a,4)*((-3. + 4.*r[jr])*pow(a,2) + (-5. + 4.*r[jr])*pow(r[jr],2))*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - 1.*r[jr]*pow(a,2)*(pow(a,4)*(5. - 15.*pow(z[jz],2) + 6.*pow(z[jz],4)) - 1.*pow(r[jr],3)*(-30. + 27.*r[jr] + 3.*(20. - 19.*r[jr])*pow(z[jz],2) + 2.*(1. + r[jr])*pow(z[jz],4)) + 2.*r[jr]*pow(a,2)*(7. - 11.*r[jr] + 3.*(-4. + 7.*r[jr])*pow(z[jz],2) + (-3. + 2.*r[jr])*pow(z[jz],4))) - I*z[jz]*pow(a,3)*(pow(a,4)*(-1. + 3.*pow(z[jz],2) - 6.*pow(z[jz],4)) + 2.*r[jr]*pow(a,2)*(-5. + 7.*r[jr] + (21. - 29.*r[jr])*pow(z[jz],2) + 6.*(-2. + 3.*r[jr])*pow(z[jz],4)) + pow(r[jr],3)*(3.*(-6. + 5.*r[jr]) + (74. - 61.*r[jr])*pow(z[jz],2) + 6.*(-8. + 7.*r[jr])*pow(z[jz],4))))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(pow(a,2) + pow(r[jr],2))*(8.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(r[jr],4)*(1. + 3.*pow(z[jz],2)) + pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*pow(a,2)*pow(r[jr],2)*(4. - 9.*pow(z[jz],2) + pow(z[jz],4)) + 4.*I*r[jr]*z[jz]*pow(a,3)*(1. - 4.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2))*(8.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(r[jr],4)*(1. + 3.*pow(z[jz],2)) + pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*pow(a,2)*pow(r[jr],2)*(4. - 9.*pow(z[jz],2) + pow(z[jz],4)) + 4.*I*r[jr]*z[jz]*pow(a,3)*(1. - 4.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*(8.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(r[jr],4)*(1. + 3.*pow(z[jz],2)) + pow(a,4)*(1. + 3.*pow(z[jz],2))*pow(z[jz],4) - 2.*pow(a,2)*pow(r[jr],2)*(4. - 9.*pow(z[jz],2) + pow(z[jz],4)) + 4.*I*r[jr]*z[jz]*pow(a,3)*(1. - 4.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(pow(a,2) + pow(r[jr],2),2)*pow(-1. + pow(z[jz],2),-3);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_22(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (I*r[jr]*z[jz] + a*(-1. + 2.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*(r[jr]*z[jz] + I*(a - 2.*a*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (I*r[jr] + a*z[jz])*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.5*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*pow(a,2)*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*a*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*a*(I*r[jr] + a*z[jz])*pow(r[jr] - I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.5*pow(a,2)*pow(r[jr] + I*a*z[jz],-2);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_22_dz(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (Complex(0.,-2.)*r[jr]*pow(a,2)*pow(z[jz],2)*(-1. + 3.*pow(z[jz],2)) - 4.*a*pow(r[jr],2)*pow(z[jz],3) - I*r[jr]*(pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) + pow(a,2)*(1. - 3.*pow(z[jz],2) + 2.*pow(z[jz],4))) - 1.*a*z[jz]*(-1.*pow(r[jr],2)*(-1. + pow(z[jz],2)) + pow(a,2)*(3. - 9.*pow(z[jz],2) + 10.*pow(z[jz],4))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-2.*r[jr]*z[jz] + I*(a - 3.*a*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = ((-I)*a*z[jz]*pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) + pow(r[jr],3)*(1. + 3.*pow(z[jz],2)) - 1.*r[jr]*pow(a,2)*(-1. + 5.*pow(z[jz],2) - 8.*pow(z[jz],4)) - I*z[jz]*pow(a,3)*(3. - 9.*pow(z[jz],2) + 10.*pow(z[jz],4)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (Complex(0.,-4.)*z[jz]*pow(r[jr],3) + I*r[jr]*z[jz]*pow(a,2)*(4. - 8.*pow(z[jz],2)) + a*pow(r[jr],2)*(-2. - 2.*pow(z[jz],2)) - 2.*pow(a,3)*pow(z[jz],2)*(-1. + 3.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (2.*r[jr]*z[jz] + I*a*(-1. + 3.*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*(r[jr] - Complex(0.,3.)*a*z[jz])*pow(a,3)*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*a*(r[jr]*z[jz] + I*a*(-1. + 2.*pow(z[jz],2)))*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(2.*a*pow(r[jr],2) + 2.*I*z[jz]*pow(r[jr],3) - 2.*pow(a,3)*(1. - 2.*pow(z[jz],2))*pow(z[jz],2) + I*r[jr]*z[jz]*pow(a,2)*(-4. + 6.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-I)*pow(a,3)*pow(r[jr] + I*a*z[jz],-3);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_22_dz2(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(6.*I*z[jz]*pow(r[jr],5)*(1. + pow(z[jz],2)) + 6.*a*pow(r[jr],4)*pow(z[jz],2)*(1. + pow(z[jz],2)) + 4.*I*pow(a,2)*pow(r[jr],3)*(-1. + pow(z[jz],2))*pow(z[jz],3) + Complex(0.,3.)*r[jr]*pow(a,4)*pow(z[jz],3)*(1. - 4.*pow(z[jz],2) + 7.*pow(z[jz],4)) - I*r[jr]*z[jz]*pow(a,2)*(pow(r[jr],2)*(5. - 10.*pow(z[jz],2) - 19.*pow(z[jz],4)) + 4.*pow(a,2)*(1. - 4.*pow(z[jz],2) + 6.*pow(z[jz],4) - 3.*pow(z[jz],6))) + pow(a,3)*pow(z[jz],2)*(pow(r[jr],2)*(-3. + 10.*pow(z[jz],2) - 7.*pow(z[jz],4)) + 6.*pow(a,2)*(-1. + 4.*pow(z[jz],2) - 6.*pow(z[jz],4) + 5.*pow(z[jz],6))) + a*pow(r[jr],2)*(pow(r[jr],2)*(1. + 2.*pow(z[jz],2) - 3.*pow(z[jz],4)) + 2.*pow(a,2)*(1. - 4.*pow(z[jz],2) + 5.*pow(z[jz],4) + 10.*pow(z[jz],6))))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (2.*pow(r[jr],2)*(1. + 5.*pow(z[jz],2)) + 4.*I*a*r[jr]*z[jz]*(-1. + 7.*pow(z[jz],2)) - 3.*pow(a,2)*(1. - 4.*pow(z[jz],2) + 7.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(-6.*z[jz]*pow(r[jr],5)*(1. + pow(z[jz],2)) + I*a*pow(r[jr],4)*(1. + 8.*pow(z[jz],2) + 3.*pow(z[jz],4)) + I*pow(a,3)*pow(r[jr],2)*(2. + pow(z[jz],2))*(1. - 6.*pow(z[jz],2) + 13.*pow(z[jz],4)) - 1.*z[jz]*pow(a,2)*pow(r[jr],3)*(-5. + 6.*pow(z[jz],2) + 23.*pow(z[jz],4)) + r[jr]*z[jz]*pow(a,4)*(4. - 19.*pow(z[jz],2) + 36.*pow(z[jz],4) - 33.*pow(z[jz],6)) + 6.*I*pow(a,5)*pow(z[jz],2)*(-1. + 4.*pow(z[jz],2) - 6.*pow(z[jz],4) + 5.*pow(z[jz],6)))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(6.*a*z[jz]*pow(r[jr],4)*(1. + pow(z[jz],2)) - 1.*a*z[jz]*pow(r[jr],2)*(4.*pow(r[jr],2)*(-1. + pow(z[jz],2)) + pow(a,2)*(5. - 10.*pow(z[jz],2) - 19.*pow(z[jz],4))) - I*(pow(a,2)*pow(r[jr],3)*(2. + 4.*pow(z[jz],2) - 30.*pow(z[jz],4)) + r[jr]*pow(a,4)*pow(z[jz],2)*(-3. + 10.*pow(z[jz],2) - 7.*pow(z[jz],4))) - 4.*z[jz]*pow(a,3)*pow(r[jr],2)*(1. - 3.*pow(z[jz],2) + 2.*pow(z[jz],4)) + 2.*I*r[jr]*pow(a,4)*pow(z[jz],2)*(3. - 11.*pow(z[jz],2) + 14.*pow(z[jz],4)) + I*(2.*pow(r[jr],5)*(1. + 5.*pow(z[jz],2)) + pow(a,2)*pow(r[jr],3)*(-1. - 2.*pow(z[jz],2) + 3.*pow(z[jz],4))) + 3.*pow(a,5)*(pow(z[jz],3) - 4.*pow(z[jz],5) + 7.*pow(z[jz],7)))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (4.*I*a*r[jr]*z[jz]*(1. - 7.*pow(z[jz],2)) - 2.*pow(r[jr],2)*(1. + 5.*pow(z[jz],2)) + 3.*pow(a,2)*(1. - 4.*pow(z[jz],2) + 7.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-4);

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 4.*pow(a,4)*(Complex(0.,-2.)*a*r[jr]*z[jz] + pow(r[jr],2) - 3.*pow(a,2)*pow(z[jz],2))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*a*(pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) + 2.*I*a*r[jr]*z[jz]*(-1. + 5.*pow(z[jz],2)) - 1.*pow(a,2)*(3. - 9.*pow(z[jz],2) + 10.*pow(z[jz],4)))*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*a*(a*z[jz]*pow(r[jr],4)*(3. + pow(z[jz],2)) - I*r[jr]*pow(a,2)*(pow(r[jr],2)*(2. - 10.*pow(z[jz],4)) + pow(a,2)*pow(z[jz],2)*(-3. + 8.*pow(z[jz],2) - 5.*pow(z[jz],4))) - 1.*a*z[jz]*pow(r[jr],2)*(2.*pow(r[jr],2)*(-1. + pow(z[jz],2)) + pow(a,2)*(5. - 10.*pow(z[jz],2) - 3.*pow(z[jz],4))) + I*pow(r[jr],3)*(pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) + pow(a,2)*(-1. + pow(z[jz],4))) + I*r[jr]*pow(a,4)*pow(z[jz],2)*(6. - 17.*pow(z[jz],2) + 15.*pow(z[jz],4)) + z[jz]*pow(a,3)*(-2.*pow(r[jr],2)*(2. - 5.*pow(z[jz],2) + 3.*pow(z[jz],4)) + pow(a,2)*pow(z[jz],2)*(3. - 9.*pow(z[jz],2) + 10.*pow(z[jz],4))))*pow(r[jr] - I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -3.*pow(a,4)*pow(r[jr] + I*a*z[jz],-4);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_24(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.4142135623730951*(I*r[jr]*((-2. + r[jr])*r[jr] + pow(a,2)) + z[jz]*pow(a,3)*(-1. + pow(z[jz],2)))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,0.7071067811865475)*a*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*r[jr]*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*z[jz]*pow(a,2)*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-0.5);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = Complex(0.,0.7071067811865475)*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -0.7071067811865475*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*a*((-I)*r[jr]*((-2. + r[jr])*r[jr] + pow(a,2)) + a*z[jz]*(pow(a,2) + pow(r[jr],2)))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-0.5);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*(-1.*pow(r[jr],2) + pow(a,2)*(-2. + pow(z[jz],2)))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (pow(a,2) + pow(r[jr],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*a*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow(2. - 2.*pow(z[jz],2),-0.5);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*a*(pow(a,2) + pow(r[jr],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow(r[jr] + I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_24_dz(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-2.5)*(pow(a,3)*pow(r[jr],2)*(-1. + pow(z[jz],2)) + Complex(0.,3.)*r[jr]*z[jz]*pow(a,2)*((-2. + r[jr])*r[jr] + pow(a,2))*(-1. + 2.*pow(z[jz],2)) - 1.*a*(-1. + pow(z[jz],2))*(pow(a,2)*pow(r[jr],2) + (-2. + r[jr])*pow(r[jr],3) + pow(a,4)*(2. - 3.*pow(z[jz],2))*pow(z[jz],2)) + I*r[jr]*z[jz]*(3.*pow(a,2)*pow(r[jr],2) + 3.*(-2. + r[jr])*pow(r[jr],3) + pow(a,4)*pow(-1. + pow(z[jz],2),2)));

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*a*((-I)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(a,2)*pow(z[jz],2)*(-1. + 4.*pow(z[jz],2)) + pow(r[jr],2)*(-3. + 6.*pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*r[jr]*(-3.*z[jz]*pow(r[jr],2) - I*a*r[jr]*(-1. + pow(z[jz],2)) - 3.*z[jz]*pow(a,2)*(-1. + 2.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(1. - 1.*pow(z[jz],2),-2.5);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(I*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2)*(3. - 6.*pow(z[jz],2)) + pow(a,2)*(1. - 4.*pow(z[jz],2))*pow(z[jz],2)) + z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.4142135623730951*pow(a,2)*(pow(r[jr],2) + I*a*r[jr]*z[jz]*(-1. + pow(z[jz],2)) - 1.*pow(a,2)*(2. - 3.*pow(z[jz],2))*pow(z[jz],2))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*((-I)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(a,2)*pow(z[jz],2)*(-1. + 4.*pow(z[jz],2)) + pow(r[jr],2)*(-3. + 6.*pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (I*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(r[jr],2)*(3. - 6.*pow(z[jz],2)) + pow(a,2)*(1. - 4.*pow(z[jz],2))*pow(z[jz],2)) + z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.4142135623730951*a*((-I)*(-2. + r[jr])*z[jz]*pow(r[jr],4) + a*pow(r[jr],3)*(r[jr] + (-2. + r[jr])*(-1. + pow(z[jz],2))) + I*r[jr]*z[jz]*pow(a,4)*(2. - 3.*pow(z[jz],2)) - 1.*pow(a,5)*(2. - 3.*pow(z[jz],2))*pow(z[jz],2) + I*z[jz]*pow(a,2)*pow(r[jr],2)*(r[jr]*(-2. + pow(z[jz],2)) - 1.*(-2. + r[jr])*(-3. + 4.*pow(z[jz],2))) + pow(a,3)*pow(r[jr],2)*(-1.*pow(z[jz],2) + 3.*pow(z[jz],4)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-1.5);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*(I*a*r[jr]*(-1.*pow(r[jr],2) + pow(a,2)*(-2. + pow(z[jz],2)))*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(6.*pow(a,2)*pow(r[jr],2)*(-1. + 2.*pow(z[jz],2)) + 3.*pow(r[jr],4)*(-1. + 2.*pow(z[jz],2)) + pow(a,4)*pow(z[jz],2)*(-2. + 7.*pow(z[jz],2) - 2.*pow(z[jz],4))) - 1.*z[jz]*pow(r[jr],2)*(6.*pow(a,2)*pow(r[jr],2) + 3.*pow(r[jr],4) + pow(a,4)*(4. - 2.*pow(z[jz],2) + pow(z[jz],4))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.*(pow(a,2) + pow(r[jr],2))*((-I)*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(r[jr],2)*(3.*pow(r[jr],2) + pow(a,2)*(2. + pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*(pow(a,2)*pow(z[jz],2)*(-1. + 4.*pow(z[jz],2)) + pow(r[jr],2)*(-3. + 6.*pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*a*(z[jz]*pow(r[jr],2)*(pow(r[jr],2) - 1.*pow(a,2)*(-2. + pow(z[jz],2))) - 1.*z[jz]*pow(a,2)*pow(r[jr],2)*(3. - 4.*pow(z[jz],2)) + I*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*pow(a,4)*(pow(z[jz],3) - 2.*pow(z[jz],5)))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-I)*a*(pow(a,2) + pow(r[jr],2))*(-1.*z[jz]*pow(r[jr],4) + z[jz]*pow(a,2)*pow(r[jr],2)*(-2. + pow(z[jz],2)) - I*a*r[jr]*(-1. + pow(z[jz],2))*(pow(r[jr],2) + pow(a,2)*pow(z[jz],2)) - 1.*z[jz]*pow(a,2)*(pow(a,2)*pow(z[jz],2)*(-1. + 2.*pow(z[jz],2)) + pow(r[jr],2)*(-3. + 4.*pow(z[jz],2))))*pow(r[jr] - I*a*z[jz],-2)*pow(r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-1);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_24_dz2(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -1.4142135623730951*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-3.5)*(3.*z[jz]*pow(a,3)*pow(r[jr],4)*(-1. + pow(z[jz],2)) - 2.*a*z[jz]*pow(r[jr],2)*(-1. + pow(z[jz],2))*(3.*pow(a,2)*pow(r[jr],2) + 3.*(-2. + r[jr])*pow(r[jr],3) + pow(a,4)*(5. + 2.*(-5. + pow(z[jz],2))*pow(z[jz],2))) + Complex(0.,3.)*r[jr]*((-2. + r[jr])*r[jr] + pow(a,2))*pow(a,4)*pow(z[jz],2)*(4. - 13.*pow(z[jz],2) + 14.*pow(z[jz],4)) + z[jz]*pow(a,3)*(-1. + pow(z[jz],2))*(2.*pow(a,2)*pow(r[jr],2)*(4. - 7.*pow(z[jz],2)) - 2.*(-2. + r[jr])*pow(r[jr],3)*(-4. + 7.*pow(z[jz],2)) + 3.*pow(a,4)*pow(z[jz],2)*(2. - 5.*pow(z[jz],2) + 4.*pow(z[jz],4))) + I*pow(r[jr],3)*(3.*pow(a,2)*pow(r[jr],2)*(1. + 4.*pow(z[jz],2)) + 3.*(-2. + r[jr])*pow(r[jr],3)*(1. + 4.*pow(z[jz],2)) + 2.*pow(a,4)*pow(-1. + pow(z[jz],2),2)) - 2.*I*r[jr]*pow(a,2)*(pow(a,2)*pow(r[jr],2)*(2. + 2.*pow(z[jz],2) - 19.*pow(z[jz],4)) - 1.*(-2. + r[jr])*pow(r[jr],3)*(-2. - 2.*pow(z[jz],2) + 19.*pow(z[jz],4)) - 1.*pow(a,4)*pow(z[jz],2)*(-3. + 4.*pow(z[jz],2))*pow(-1. + pow(z[jz],2),2)));

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-I)*a*(6.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) - 2.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(2. - 8.*pow(z[jz],2))*(-1. + pow(z[jz],2)) + 3.*pow(r[jr],6)*(1. + 4.*pow(z[jz],2)) + 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 7.*pow(z[jz],2) + 5.*pow(z[jz],4)) + pow(a,6)*pow(z[jz],4)*(2. - 7.*pow(z[jz],2) + 20.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],4)*(-2. + 7.*pow(z[jz],2) + 40.*pow(z[jz],4)) - 3.*pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(-4. + 13.*pow(z[jz],2) - 14.*pow(z[jz],4) - 2.*(-2. + 6.*pow(z[jz],2) + pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*r[jr]*(6.*I*a*z[jz]*pow(r[jr],3)*(-1. + pow(z[jz],2)) + 3.*pow(r[jr],4)*(1. + 4.*pow(z[jz],2)) - 2.*pow(a,2)*pow(r[jr],2)*(2. + 2.*pow(z[jz],2) - 19.*pow(z[jz],4)) + 2.*I*r[jr]*z[jz]*pow(a,3)*(4. - 11.*pow(z[jz],2) + 7.*pow(z[jz],4)) + 3.*pow(a,4)*pow(z[jz],2)*(4. - 13.*pow(z[jz],2) + 14.*pow(z[jz],4)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-3.5);

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = a*(6.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) - 2.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(2. - 8.*pow(z[jz],2))*(-1. + pow(z[jz],2)) + 3.*pow(r[jr],6)*(1. + 4.*pow(z[jz],2)) + 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 7.*pow(z[jz],2) + 5.*pow(z[jz],4)) + pow(a,6)*pow(z[jz],4)*(2. - 7.*pow(z[jz],2) + 20.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],4)*(-2. + 7.*pow(z[jz],2) + 40.*pow(z[jz],4)) - 3.*pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(-4. + 13.*pow(z[jz],2) - 14.*pow(z[jz],4) - 2.*(-2. + 6.*pow(z[jz],2) + pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*pow(a,2)*(3.*z[jz]*pow(r[jr],4) + 2.*I*a*pow(r[jr],3)*(-1. + pow(z[jz],2)) - 2.*z[jz]*pow(a,2)*pow(r[jr],2)*(5. + 2.*(-5. + pow(z[jz],2))*pow(z[jz],2)) + 2.*I*r[jr]*pow(a,3)*pow(z[jz],2)*(3. - 7.*pow(z[jz],2) + 4.*pow(z[jz],4)) + 3.*pow(a,4)*pow(z[jz],3)*(2. - 5.*pow(z[jz],2) + 4.*pow(z[jz],4)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(1. - 1.*pow(z[jz],2),-2.5);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (-I)*(6.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) - 2.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(2. - 8.*pow(z[jz],2))*(-1. + pow(z[jz],2)) + 3.*pow(r[jr],6)*(1. + 4.*pow(z[jz],2)) + 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 7.*pow(z[jz],2) + 5.*pow(z[jz],4)) + pow(a,6)*pow(z[jz],4)*(2. - 7.*pow(z[jz],2) + 20.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],4)*(-2. + 7.*pow(z[jz],2) + 40.*pow(z[jz],4)) - 3.*pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(-4. + 13.*pow(z[jz],2) - 14.*pow(z[jz],4) - 2.*(-2. + 6.*pow(z[jz],2) + pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (6.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) - 2.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(2. - 8.*pow(z[jz],2))*(-1. + pow(z[jz],2)) + 3.*pow(r[jr],6)*(1. + 4.*pow(z[jz],2)) + 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 7.*pow(z[jz],2) + 5.*pow(z[jz],4)) + pow(a,6)*pow(z[jz],4)*(2. - 7.*pow(z[jz],2) + 20.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],4)*(-2. + 7.*pow(z[jz],2) + 40.*pow(z[jz],4)) - 3.*pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(-4. + 13.*pow(z[jz],2) - 14.*pow(z[jz],4) - 2.*(-2. + 6.*pow(z[jz],2) + pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 1.4142135623730951*a*(3.*a*z[jz]*(pow(a,2) + pow(r[jr],2))*pow(r[jr],4) - I*pow(r[jr],3)*(3.*pow(a,2)*pow(r[jr],2) - 2.*pow(a,4)*(-1. + pow(z[jz],2)) + (-2. + r[jr])*pow(r[jr],3)*(1. + 2.*pow(z[jz],2))) - I*r[jr]*((-2. + r[jr])*r[jr] + pow(a,2))*pow(a,4)*pow(z[jz],2)*(12. - 29.*pow(z[jz],2) + 20.*pow(z[jz],4)) - 2.*a*z[jz]*pow(r[jr],2)*(-1.*(-2. + r[jr])*pow(r[jr],3)*(-1. + pow(z[jz],2)) + pow(a,4)*(5. + 2.*(-5. + pow(z[jz],2))*pow(z[jz],2)) + pow(a,2)*pow(r[jr],2)*(6. - 11.*pow(z[jz],2) + 2.*pow(z[jz],4))) - 2.*I*pow(a,2)*(r[jr]*pow(a,4)*pow(z[jz],2)*(-3. + 7.*pow(z[jz],2) - 4.*pow(z[jz],4)) + (-2. + r[jr])*pow(r[jr],4)*(-2. + 2.*pow(z[jz],2) + 3.*pow(z[jz],4)) - 1.*pow(a,2)*pow(r[jr],3)*(2. + pow(z[jz],2) - 10.*pow(z[jz],4) + 4.*pow(z[jz],6))) + z[jz]*pow(a,3)*(3.*pow(a,4)*pow(z[jz],2)*(2. - 5.*pow(z[jz],2) + 4.*pow(z[jz],4)) + 2.*(-2. + r[jr])*pow(r[jr],3)*(4. - 9.*pow(z[jz],2) + 5.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],2)*(8. - 12.*pow(z[jz],2) - 5.*pow(z[jz],4) + 12.*pow(z[jz],6))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(1. - 1.*pow(z[jz],2),-2.5);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*((-6.*I)*a*z[jz]*pow(r[jr],7)*(-1. + pow(z[jz],2)) + 2.*I*z[jz]*pow(a,3)*pow(r[jr],5)*(-2. - 7.*pow(z[jz],2))*(-1. + pow(z[jz],2)) - 3.*pow(r[jr],8)*(1. + 4.*pow(z[jz],2)) + 2.*pow(a,2)*pow(r[jr],6)*(-1. - 10.*pow(z[jz],2) - 19.*pow(z[jz],4)) - 2.*I*z[jz]*pow(a,5)*pow(r[jr],3)*(-1. + pow(z[jz],2))*(-4. + 12.*pow(z[jz],2) + pow(z[jz],4)) - 1.*pow(a,4)*pow(r[jr],4)*(-4. + 2.*pow(z[jz],2) + 77.*pow(z[jz],4) + 3.*pow(z[jz],2)*(4. - 13.*pow(z[jz],2) + 14.*pow(z[jz],4))) + 2.*pow(a,6)*pow(r[jr],2)*pow(z[jz],2)*(12. - 38.*pow(z[jz],2) + 16.*pow(z[jz],4) - 3.*(4. - 13.*pow(z[jz],2) + 14.*pow(z[jz],4)) - 5.*pow(z[jz],6)) - 2.*I*r[jr]*pow(a,7)*pow(z[jz],3)*(4. - 14.*pow(z[jz],2) + 13.*pow(z[jz],4) - 3.*pow(z[jz],6)) + pow(a,8)*pow(z[jz],4)*(-4. + 14.*pow(z[jz],2) - 31.*pow(z[jz],4) + 6.*pow(z[jz],6)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (pow(a,2) + pow(r[jr],2))*(6.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) - 2.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(2. - 8.*pow(z[jz],2))*(-1. + pow(z[jz],2)) + 3.*pow(r[jr],6)*(1. + 4.*pow(z[jz],2)) + 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 7.*pow(z[jz],2) + 5.*pow(z[jz],4)) + pow(a,6)*pow(z[jz],4)*(2. - 7.*pow(z[jz],2) + 20.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],4)*(-2. + 7.*pow(z[jz],2) + 40.*pow(z[jz],4)) - 3.*pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(-4. + 13.*pow(z[jz],2) - 14.*pow(z[jz],4) - 2.*(-2. + 6.*pow(z[jz],2) + pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*a*(2.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) - 2.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(2. - 4.*pow(z[jz],2))*(-1. + pow(z[jz],2)) + pow(r[jr],6)*(1. + 2.*pow(z[jz],2)) - 1.*pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(pow(z[jz],2) - 10.*pow(z[jz],4)) + 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 5.*pow(z[jz],2) + 3.*pow(z[jz],4)) + pow(a,6)*pow(z[jz],4)*(2. - 5.*pow(z[jz],2) + 6.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],4)*(-2. + 5.*pow(z[jz],2) + 6.*pow(z[jz],4)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = I*a*(pow(a,2) + pow(r[jr],2))*(2.*I*a*z[jz]*pow(r[jr],5)*(-1. + pow(z[jz],2)) - 2.*I*z[jz]*pow(a,3)*pow(r[jr],3)*(2. - 4.*pow(z[jz],2))*(-1. + pow(z[jz],2)) + pow(r[jr],6)*(1. + 2.*pow(z[jz],2)) - 1.*pow(a,4)*pow(r[jr],2)*pow(z[jz],2)*(pow(z[jz],2) - 10.*pow(z[jz],4)) + 2.*I*r[jr]*pow(a,5)*pow(z[jz],3)*(2. - 5.*pow(z[jz],2) + 3.*pow(z[jz],4)) + pow(a,6)*pow(z[jz],4)*(2. - 5.*pow(z[jz],2) + 6.*pow(z[jz],4)) + pow(a,2)*pow(r[jr],4)*(-2. + 5.*pow(z[jz],2) + 6.*pow(z[jz],4)))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(r[jr] + I*a*z[jz],-4)*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(2. - 2.*pow(z[jz],2),-0.5)*pow(-1. + pow(z[jz],2),-2);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_44(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*a*(r[jr]*(-3. + 2.*r[jr]) - I*a*(-1. + r[jr])*z[jz] + pow(a,2))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = pow(a,2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*pow(r[jr] - I*a*z[jz],-1)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*a*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(r[jr]*(-3. + 2.*r[jr])*pow(a,2) + I*z[jz]*pow(a,3) + pow(a,4) - I*a*z[jz]*pow(r[jr],2) + (-1. + r[jr])*pow(r[jr],3))*pow(-1.*r[jr] + I*a*z[jz],-1)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*a*(pow(a,2) + pow(r[jr],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(pow(a,2) + pow(r[jr],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(-1. + pow(z[jz],2),-1);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(pow(a,2) + pow(r[jr],2),2)*pow(-1. + pow(z[jz],2),-1);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_44_dz(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*a*(2.*r[jr]*z[jz]*(r[jr]*(-3. + 2.*r[jr]) + pow(a,2)) + I*a*((-2. + r[jr])*r[jr] + pow(a,2) + ((10. - 7.*r[jr])*r[jr] - 3.*pow(a,2))*pow(z[jz],2)) - 2.*(-1. + r[jr])*pow(a,2)*pow(z[jz],3))*pow(r[jr] - I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*z[jz]*pow(a,2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (4.*r[jr]*z[jz] + I*a*(2. - 6.*pow(z[jz],2)))*pow(r[jr] - I*a*z[jz],-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -4.*a*z[jz]*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*z[jz]*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (2.*I*a*((-2. + r[jr])*r[jr] + pow(a,2))*(pow(a,2) + pow(r[jr],2)) + 4.*r[jr]*z[jz]*(r[jr]*(-3. + 2.*r[jr])*pow(a,2) + pow(a,4) + (-1. + r[jr])*pow(r[jr],3)) - 2.*I*a*(2.*r[jr]*(-5. + 3.*r[jr])*pow(a,2) + 3.*pow(a,4) + (-2. + 3.*r[jr])*pow(r[jr],3))*pow(z[jz],2) + 4.*(a - 1.*r[jr])*(a + r[jr])*pow(a,2)*pow(z[jz],3))*pow(r[jr] - I*a*z[jz],-2)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -4.*a*z[jz]*(pow(a,2) + pow(r[jr],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -4.*z[jz]*(pow(a,2) + pow(r[jr],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(-1. + pow(z[jz],2),-2);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -2.*z[jz]*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(pow(a,2) + pow(r[jr],2),2)*pow(-1. + pow(z[jz],2),-2);
// 		}
// 	}

// 	return coeffs;
// }

// ComplexTensor metric_coefficients_IRG_44_dz2(double a, Vector r, Vector z){
// 	int componentNum = 3*3*3*3;
// 	int ai, bi, ci, di, pos;
// 	ComplexTensor coeffs(componentNum, ComplexMatrix(r.size(), ComplexVector(z.size())));

// 	for(size_t jr = 0; jr < r.size(); jr++){
// 		for(size_t jz = 0; jz < z.size(); jz++){
// 		ai = 0, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 4.*a*((r[jr]*(-3. + 2.*r[jr]) + pow(a,2))*pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) + I*a*r[jr]*z[jz]*(-3.*(-1. + r[jr])*r[jr] + ((25. - 17.*r[jr])*r[jr] - 8.*pow(a,2))*pow(z[jz],2)) + I*(-1. + r[jr])*pow(a,3)*(1. + 3.*pow(z[jz],2))*pow(z[jz],3) - 1.*pow(a,2)*(pow(a,2)*(1. - 3.*pow(z[jz],2) + 6.*pow(z[jz],4)) + r[jr]*(-2. + r[jr] + 3.*pow(z[jz],2) + 3.*(-7. + 5.*r[jr])*pow(z[jz],4))))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 0, di = 2;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*pow(a,2)*(1. + 3.*pow(z[jz],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 1, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 0, ci = 2, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = -4.*(pow(r[jr],2)*(1. + 3.*pow(z[jz],2)) - 8.*I*a*r[jr]*pow(z[jz],3) - 1.*pow(a,2)*(1. - 3.*pow(z[jz],2) + 6.*pow(z[jz],4)))*pow(r[jr] - I*a*z[jz],-3)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 4.*a*(1. + 3.*pow(z[jz],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 0, bi = 1, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 0, bi = 2, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = (2. + 6.*pow(z[jz],2))*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*(pow(I*a - 1.*r[jr],-1)*(I*a*(a - 1.*r[jr])*(a + r[jr]) + r[jr]*(-3. + 2.*r[jr])*pow(a,2) + pow(a,4) + (-1. + r[jr])*pow(r[jr],3))*pow(-1. + z[jz],-3) + (r[jr]*(-3. + 2.*r[jr])*pow(a,2) + pow(a,4) + I*a*(-1.*pow(a,2) + pow(r[jr],2)) + (-1. + r[jr])*pow(r[jr],3))*pow(I*a + r[jr],-1)*pow(1. + z[jz],-3) - 2.*pow(a,4)*(pow(a,2) + pow(r[jr],2))*(-2.*r[jr] + pow(a,2) + pow(r[jr],2))*pow(-1.*r[jr] + I*a*z[jz],-3)*pow(-1.*pow(a,2) - 1.*pow(r[jr],2),-1));

// 		ai = 1, bi = 0, ci = 0, di = 1;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 4.*a*(pow(a,2) + pow(r[jr],2))*(1. + 3.*pow(z[jz],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 1, bi = 0, ci = 1, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 0.;

// 		ai = 1, bi = 1, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 4.*(pow(a,2) + pow(r[jr],2))*(1. + 3.*pow(z[jz],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-1)*pow(-1. + pow(z[jz],2),-3);

// 		ai = 2, bi = 0, ci = 0, di = 0;
// 		pos = 3*(3*((3*ai) + bi) + ci) + di;
// 		coeffs[pos][jr][jz] = 2.*(1. + 3.*pow(z[jz],2))*pow((-2. + r[jr])*r[jr] + pow(a,2),-2)*pow(pow(a,2) + pow(r[jr],2),2)*pow(-1. + pow(z[jz],2),-3);
// 		}
// 	}

// 	return coeffs;
// }

Complex u1_subfunction(double a, double En, double Lz, double, double rp, double zp, double deltaUr){
  return 0.5*pow(rp*rp + a*a*zp*zp, -1)*(En*(rp*rp + a*a) - a*Lz + deltaUr);
}

Complex u1_dz_subfunction(double a, double En, double Lz, double, double rp, double zp, double deltaUr){
  return -(a*a*zp)*pow(rp*rp + a*a*zp*zp, -2)*(En*(rp*rp + a*a) - a*Lz + deltaUr);
}

Complex u1_dz2_subfunction(double a, double En, double Lz, double, double rp, double zp, double deltaUr){
  return -a*a*(rp*rp - 3.*a*a*zp*zp)*pow(rp*rp + a*a*zp*zp, -3)*(En*(rp*rp + a*a) - a*Lz + deltaUr);
}

Complex u2_subfunction(double a, double En, double Lz, double, double rp, double, double deltaUr){
  return pow(rp*rp - 2.*rp + a*a, -1)*(En*(rp*rp + a*a) - a*Lz - deltaUr);
}

Complex u2_dz_subfunction(double, double, double, double, double, double, double){
  return 0.;
}

Complex u2_dz2_subfunction(double, double, double, double, double, double, double){
  return 0.;
}

// Complex u3_subfunction(double a, double En, double Lz, double, double r, double z, double Uz){
//   return (-I)*Uz*pow(I*r + a*z,-1)*pow(2. - 2.*pow(z,2),-0.5) + (Lz + a*En*(-1. + pow(z,2)))*pow(I*r + a*z,-1)*pow(2. - 2.*pow(z,2),-0.5);
// }

// Complex u3_dz_subfunction(double a, double En, double Lz, double, double r, double z, double Uz, double Uzdz){
//   return (Uz*(r*z + I*a*(1. - 2.*pow(z,2)))*pow(I*r + a*z,-2)*pow(1. - 1.*pow(z,2),-1.5) - 1.*((-I)*Lz*r*z + En*pow(a,2)*(-1. + pow(z,2)) + a*(Lz + I*En*r*z*(-1. + pow(z,2)) - 2.*Lz*pow(z,2)))*pow(I*r + a*z,-2)*pow(1. - 1.*pow(z,2),-1.5) - 1.*Uzdz*pow(r - I*a*z,-1)*pow(1. - 1.*pow(z,2),-0.5))/sqrt(2.);
// }

// Complex u3_dz2_subfunction(double a, double En, double Lz, double, double r, double z, double Uz, double Uzdz, double Uzdz2){
//   return (Uz*(I*pow(r,2)*(1. + 2.*pow(z,2)) + 6.*a*r*pow(z,3) - I*pow(a,2)*(2. - 5.*pow(z,2) + 6.*pow(z,4)))*pow(I*r + a*z,-3)*pow(1. - 1.*pow(z,2),-2.5) + (-1.*Lz*pow(r,2)*(1. + 2.*pow(z,2)) + a*r*(En*r*(-1. + pow(z,2)) + 6.*I*Lz*pow(z,3)) + En*pow(a,3)*(-2. + 5.*pow(z,2) - 3.*pow(z,4)) + pow(a,2)*((-2.*I)*En*r*(-1. + pow(z,2))*pow(z,3) + Lz*(2. - 5.*pow(z,2) + 6.*pow(z,4))))*pow(I*r + a*z,-3)*pow(1. - 1.*pow(z,2),-2.5) + 2.*Uzdz*(r*z + I*a*(1. - 2.*pow(z,2)))*pow(I*r + a*z,-2)*pow(1. - 1.*pow(z,2),-1.5) - 1.*Uzdz2*pow(r - I*a*z,-1)*pow(1. - 1.*pow(z,2),-0.5))/sqrt(2);
// }

Complex u3_subfunction(double a, double En, double Lz, double, double r, double z, double Uz){
  return (Lz + a*En*(-1. + pow(z,2)) + I*Uz*(-1. + pow(z,2)))*pow(I*r + a*z,-1)*pow(2. - 2.*pow(z,2),-0.5);
}

Complex u3_dz_subfunction(double a, double En, double Lz, double, double r, double z, double Uz, double Uzdz){
	return (-I*Lz*r*z + En*pow(a,2)*(-1. + pow(z,2)) + (I*Uz*(a + I*r*z) - 1.*Uzdz*(r - I*a*z)*(-1. + pow(z,2)))*(-1. + pow(z,2)) + a*(Lz + I*En*r*z*(-1. + pow(z,2)) - 2.*Lz*pow(z,2)))*pow(I*r + a*z,-2)*pow(2. - 2.*pow(z,2),-0.5)*pow(-1. + pow(z,2),-1);
}

Complex u3_dz2_subfunction(double a, double En, double Lz, double, double r, double z, double Uz, double Uzdz, double Uzdz2){
	return -I*pow(I*r + a*z,-3)*pow(2. - 2.*pow(z,2),-0.5)*pow(-1. + pow(z,2),-2)*(Uz*(-1. + pow(z,2))*(-1.*pow(r,2) + pow(a,2)*(-2. + 3.*pow(z,2)) + 2.*I*a*r*pow(z,3)) + 2.*Uzdz*(r - I*a*z)*(-I*a + r*z)*pow(-1. + pow(z,2),2) - I*(Lz*pow(r,2)*(1. + 2.*pow(z,2)) + a*r*(En*(r - 1.*r*pow(z,2)) - 6.*I*Lz*pow(z,3)) + pow(a,2)*(2.*I*En*r*(-1. + pow(z,2))*pow(z,3) + Lz*(-2. + 5.*pow(z,2) - 6.*pow(z,4))) + En*pow(a,3)*(2. - 5.*pow(z,2) + 3.*pow(z,4)) + I*Uzdz2*pow(r - I*a*z,2)*pow(-1. + pow(z,2),3)));
}

Complex u3_subfunction(double a, double En, double Lz, double Qc, double r, double z, double Uz, double Uzdz, double Uzdz2){
  return u3_subfunction(a, En, Lz, Qc, r, z, Uz);
}

Complex u3_dz_subfunction(double a, double En, double Lz, double Qc, double r, double z, double Uz, double Uzdz, double Uzdz2){
  return u3_dz_subfunction(a, En, Lz, Qc, r, z, Uz, Uzdz);
}

Complex u4_subfunction(double a, double En, double Lz, double Qc, double r, double z, double Uz){
  return std::conj(u3_subfunction(a, En, Lz, Qc, r, z, Uz));
}

Complex u4_dz_subfunction(double a, double En, double Lz, double Qc, double r, double z, double Uz, double Uzdz){
  return std::conj(u3_dz_subfunction(a, En, Lz, Qc, r, z, Uz, Uzdz));
}

Complex u4_subfunction(double a, double En, double Lz, double Qc, double r, double z, double Uz, double Uzdz, double Uzdz2){
  return u4_subfunction(a, En, Lz, Qc, r, z, Uz);
}

Complex u4_dz_subfunction(double a, double En, double Lz, double Qc, double r, double z, double Uz, double Uzdz, double Uzdz2){
  return u4_dz_subfunction(a, En, Lz, Qc, r, z, Uz, Uzdz);
}

Complex u4_dz2_subfunction(double a, double En, double Lz, double Qc, double r, double z, double Uz, double Uzdz, double Uzdz2){
  return std::conj(u3_dz2_subfunction(a, En, Lz, Qc, r, z, Uz, Uzdz, Uzdz2));
}

ComplexMatrix tetrad_velocity_radial(Complex (*subfunc)(double, double, double, double, double, double, double), double a, double En, double Lz, double Qc, Vector r, Vector z){
  int u1rSize = 2*(r.size() - 1);
  int u1zSize = 2*(z.size() - 1);
  ComplexMatrix ua(u1rSize, ComplexVector(u1zSize));
  double rp, zp;
  // first deal with double turning points
  // rmin, zmin
  rp = r[0];
  zp = z[0];
  ua[0][0] = subfunc(a, En, Lz, Qc, rp, zp, 0.);

  // rmin, zmax
  rp = r[0];
  zp = z[(z.size() - 1)];
  ua[0][(z.size() - 1)] = subfunc(a, En, Lz, Qc, rp, zp, 0.);

  // rmax, zmax
  rp = r[(r.size() - 1)];
  zp = z[(z.size() - 1)];
  ua[(r.size() - 1)][(z.size() - 1)] = subfunc(a, En, Lz, Qc, rp, zp, 0.);

  // rmax, zmin
  rp = r[(r.size() - 1)];
  zp = z[0];
  ua[(r.size() - 1)][0] = subfunc(a, En, Lz, Qc, rp, zp, 0.);

  // take into account when the radial velocity is zero but different polar velocities
  for(size_t jz = 1; jz < z.size() - 1; jz++){
    rp = r[0];
    zp = z[jz];
    ua[0][jz] = subfunc(a, En, Lz, Qc, rp, zp, 0.);
    ua[0][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, 0.);

    rp = r[(r.size() - 1)];
    zp = z[jz];
    ua[(r.size() - 1)][jz] = subfunc(a, En, Lz, Qc, rp, zp, 0.);
    ua[(r.size() - 1)][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, 0.);
  }

  // take into account when the polar velocity is zero but different radial velocities
  for(size_t jr = 1; jr < r.size() - 1; jr++){
    rp = r[jr];
    zp = z[0];
	double deltaUr = std::sqrt(std::abs(kerr_geo_Vr(a, En, Lz, Qc, rp)));
	// take into account when the radial velocity is positive
	ua[jr][0] = subfunc(a, En, Lz, Qc, rp, zp, deltaUr);
	ua[u1rSize - jr][0] = subfunc(a, En, Lz, Qc, rp, zp, -deltaUr);

	rp = r[jr];
    zp = z[z.size() - 1];
	// take into account when the radial velocity is negative
	ua[jr][z.size() - 1] = subfunc(a, En, Lz, Qc, rp, zp, deltaUr);
	ua[u1rSize - jr][z.size() - 1] = subfunc(a, En, Lz, Qc, rp, zp, -deltaUr);
  }

  // take into account different radial and polar velocities
  for(size_t jr = 1; jr < r.size() - 1; jr++){
    for(size_t jz = 1; jz < z.size() - 1; jz++){
      rp = r[jr];
      zp = z[jz];
      double deltaUr = std::sqrt(std::abs(kerr_geo_Vr(a, En, Lz, Qc, rp)));
      // take into account when the radial velocity is positive
      ua[jr][jz] = subfunc(a, En, Lz, Qc, rp, zp, deltaUr);
      ua[jr][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, deltaUr);

      // take into account when the radial velocity is negative
      ua[u1rSize - jr][jz] = subfunc(a, En, Lz, Qc, rp, zp, -deltaUr);
      ua[u1rSize - jr][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, -deltaUr);
    }
  }

  return ua;
}

ComplexMatrix tetrad_velocity_polar(Complex (*subfunc)(double, double, double, double, double, double, double, double, double), double a, double En, double Lz, double Qc, Vector r, Vector z){
  int u1rSize = 2*(r.size() - 1);
  int u1zSize = 2*(z.size() - 1);
  ComplexMatrix ua(u1rSize, ComplexVector(u1zSize));
  double rp, zp;
  // first deal with double turning points
  // rmin, zmin
  rp = r[0];
  zp = z[0];
  ua[0][0] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

  // rmin, zmax
  rp = r[0];
  zp = z[(z.size() - 1)];
  ua[0][(z.size() - 1)] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

  // rmax, zmax
  rp = r[(r.size() - 1)];
  zp = z[(z.size() - 1)];
  ua[(r.size() - 1)][(z.size() - 1)] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

  // rmax, zmin
  rp = r[(r.size() - 1)];
  zp = z[0];
  ua[(r.size() - 1)][0] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

  // take into account when the radial velocity is zero but different polar velocities
  for(size_t jz = 1; jz < z.size() - 1; jz++){
    rp = r[0];
    zp = z[jz];
	// the minus sign is because we start with utheta > 0 ==> uz < 0
	double uz = -sqrt(std::abs(kerr_geo_Vz(a, En, Lz, Qc, zp)))/(1. - zp*zp);
    ua[0][jz] = subfunc(a, En, Lz, Qc, rp, zp, uz, 0., 0.);
    ua[0][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, -uz, 0., 0.);

    rp = r[(r.size() - 1)];
    zp = z[jz];
    ua[(r.size() - 1)][jz] = subfunc(a, En, Lz, Qc, rp, zp, uz, 0., 0.);
    ua[(r.size() - 1)][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, -uz, 0., 0.);
  }

  // take into account when the polar velocity is zero but different radial velocities
  for(size_t jr = 1; jr < r.size() - 1; jr++){
    rp = r[jr];
    zp = z[0];
	// take into account when the radial velocity is positive
	ua[jr][0] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
	ua[u1rSize - jr][0] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

	rp = r[jr];
    zp = z[(z.size() - 1)];
	// take into account when the radial velocity is negative
	ua[jr][z.size() - 1] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
	ua[u1rSize - jr][z.size() - 1] = subfunc(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
  }

  // take into account different radial and polar velocities
  for(size_t jr = 1; jr < r.size() - 1; jr++){
    for(size_t jz = 1; jz < z.size() - 1; jz++){
      rp = r[jr];
      zp = z[jz];
      double uz = -sqrt(std::abs(kerr_geo_Vz(a, En, Lz, Qc, zp)))/(1. - zp*zp);
      // take into account when the radial velocity is positive
      ua[jr][jz] = subfunc(a, En, Lz, Qc, rp, zp, uz, 0., 0.);
      ua[jr][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, -uz, 0., 0.);

      // take into account when the radial velocity is negative
      ua[u1rSize - jr][jz] = subfunc(a, En, Lz, Qc, rp, zp, uz, 0., 0.);
      ua[u1rSize - jr][u1zSize - jz] = subfunc(a, En, Lz, Qc, rp, zp, -uz, 0., 0.);
    }
  }

  return ua;
}

ComplexMatrix tetrad_velocity_1(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_radial(u1_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_1_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_radial(u1_dz_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_1_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_radial(u1_dz2_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_2(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_radial(u2_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_2_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_radial(u2_dz_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_2_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_radial(u2_dz2_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_3(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_polar(u3_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_3_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_polar(u3_dz_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_3_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_polar(u3_dz2_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_4(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_polar(u4_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_4_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_polar(u4_dz_subfunction, a, En, Lz, Qc, r, z);
}

ComplexMatrix tetrad_velocity_4_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
	return tetrad_velocity_polar(u4_dz2_subfunction, a, En, Lz, Qc, r, z);
}

// ComplexMatrix tetrad_velocity_1(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u1(u1rSize, ComplexVector(u1zSize));
//   double rp, zp;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u1[0][0] = u1_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u1[(r.size() - 1)][(z.size() - 1)] = u1_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the radial velocity is zero
//   for(size_t jz = 1; jz < z.size() - 1; jz++){
//     rp = r[0];
//     zp = z[jz];
//     u1[0][jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u1[0][u1zSize - jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//     rp = r[(r.size() - 1)];
//     zp = z[jz];
//     u1[(r.size() - 1)][jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u1[(r.size() - 1)][u1zSize - jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       double deltaUr = sqrt(kerr_geo_Vr(a, En, Lz, Qc, rp));
//       // take into account when the radial velocity is positive
//       u1[jr][jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);
//       u1[jr][u1zSize - jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);

//       // take into account when the radial velocity is negative
//       u1[u1rSize - jr][jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//       u1[u1rSize - jr][u1zSize - jz] = u1_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//     }
//   }

//   return u1;
// }

// ComplexMatrix tetrad_velocity_1_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u1(u1rSize, ComplexVector(u1zSize));
//   double rp, zp;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u1[0][0] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u1[(r.size() - 1)][(z.size() - 1)] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the radial velocity is zero
//   for(size_t jz = 1; jz < z.size() - 1; jz++){
//     rp = r[0];
//     zp = z[jz];
//     u1[0][jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u1[0][u1zSize - jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//     rp = r[(r.size() - 1)];
//     zp = z[jz];
//     u1[(r.size() - 1)][jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u1[(r.size() - 1)][u1zSize - jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       double deltaUr = sqrt(kerr_geo_Vr(a, En, Lz, Qc, rp));
//       // take into account when the radial velocity is positive
//       u1[jr][jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);
//       u1[jr][u1zSize - jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);

//       // take into account when the radial velocity is negative
//       u1[u1rSize - jr][jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//       u1[u1rSize - jr][u1zSize - jz] = u1_dz_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//     }
//   }

//   return u1;
// }

// ComplexMatrix tetrad_velocity_1_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u1(u1rSize, ComplexVector(u1zSize));
//   double rp, zp;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u1[0][0] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u1[(r.size() - 1)][(z.size() - 1)] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the radial velocity is zero
//   for(size_t jz = 1; jz < z.size() - 1; jz++){
//     rp = r[0];
//     zp = z[jz];
//     u1[0][jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u1[0][u1zSize - jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//     rp = r[(r.size() - 1)];
//     zp = z[jz];
//     u1[(r.size() - 1)][jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u1[(r.size() - 1)][u1zSize - jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       double deltaUr = sqrt(kerr_geo_Vr(a, En, Lz, Qc, rp));
//       // take into account when the radial velocity is positive
//       u1[jr][jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);
//       u1[jr][u1zSize - jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);

//       // take into account when the radial velocity is negative
//       u1[u1rSize - jr][jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//       u1[u1rSize - jr][u1zSize - jz] = u1_dz2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//     }
//   }

//   return u1;
// }

// ComplexMatrix tetrad_velocity_2(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u2(u1rSize, ComplexVector(u1zSize));
//   double rp, zp;
//   // first deal with double turning points
//   // rmin, zmin
//   rp = r[0];
//   zp = z[0];
//   u2[0][0] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // rmin, zmax
//   rp = r[0];
//   zp = z[(z.size() - 1)];
//   u2[0][(z.size() - 1)] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // rmax, zmax
//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u2[(r.size() - 1)][(z.size() - 1)] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // rmax, zmin
//   rp = r[(r.size() - 1)];
//   zp = z[0];
//   u2[(r.size() - 1)][0] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the radial velocity is zero but different polar velocities
//   for(size_t jz = 1; jz < z.size() - 1; jz++){
//     rp = r[0];
//     zp = z[jz];
//     u2[0][jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u2[0][u1zSize - jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//     rp = r[(r.size() - 1)];
//     zp = z[jz];
//     u2[(r.size() - 1)][jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u2[(r.size() - 1)][u1zSize - jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   // take into account when the polar velocity is zero but different radial velocities
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     rp = r[jr];
//     zp = z[0];
// 	double deltaUr = sqrt(kerr_geo_Vr(a, En, Lz, Qc, rp));
// 	// take into account when the radial velocity is positive
// 	u2[jr][0] = u2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);
// 	u2[jr][z.size() - 1] = u2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);

// 	// take into account when the radial velocity is negative
// 	u2[u1rSize - jr][0] = u2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
// 	u2[u1rSize - jr][z.size() - 1] = u2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//   }

//   // take into account different radial and polar velocities
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       double deltaUr = sqrt(kerr_geo_Vr(a, En, Lz, Qc, rp));
//       // take into account when the radial velocity is positive
//       u2[jr][jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);
//       u2[jr][u1zSize - jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);

//       // take into account when the radial velocity is negative
//       u2[u1rSize - jr][jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//       u2[u1rSize - jr][u1zSize - jz] = u2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//     }
//   }

//   return u2;
// }

// ComplexMatrix tetrad_velocity_2_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u2(u1rSize, ComplexVector(u1zSize));
//   double rp, zp;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u2[0][0] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u2[(r.size() - 1)][(z.size() - 1)] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the radial velocity is zero
//   for(size_t jz = 1; jz < z.size() - 1; jz++){
//     rp = r[0];
//     zp = z[jz];
//     u2[0][jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u2[0][u1zSize - jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//     rp = r[(r.size() - 1)];
//     zp = z[jz];
//     u2[(r.size() - 1)][jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u2[(r.size() - 1)][u1zSize - jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       double deltaUr = sqrt(kerr_geo_Vr(a, En, Lz, Qc, rp));
//       // take into account when the radial velocity is positive
//       u2[jr][jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);
//       u2[jr][u1zSize - jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);

//       // take into account when the radial velocity is negative
//       u2[u1rSize - jr][jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//       u2[u1rSize - jr][u1zSize - jz] = u2_dz_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//     }
//   }

//   return u2;
// }

// ComplexMatrix tetrad_velocity_2_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u2(u1rSize, ComplexVector(u1zSize));
//   double rp, zp;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u2[0][0] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u2[(r.size() - 1)][(z.size() - 1)] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the radial velocity is zero
//   for(size_t jz = 1; jz < z.size() - 1; jz++){
//     rp = r[0];
//     zp = z[jz];
//     u2[0][jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u2[0][u1zSize - jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//     rp = r[(r.size() - 1)];
//     zp = z[jz];
//     u2[(r.size() - 1)][jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u2[(r.size() - 1)][u1zSize - jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       double deltaUr = sqrt(kerr_geo_Vr(a, En, Lz, Qc, rp));
//       // take into account when the radial velocity is positive
//       u2[jr][jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);
//       u2[jr][u1zSize - jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, deltaUr);

//       // take into account when the radial velocity is negative
//       u2[u1rSize - jr][jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//       u2[u1rSize - jr][u1zSize - jz] = u2_dz2_subfunction(a, En, Lz, Qc, rp, zp, -deltaUr);
//     }
//   }

//   return u2;
// }

// ComplexMatrix tetrad_velocity_3(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u3(u1rSize, ComplexVector(u1zSize));
//   double rp, zp, uz;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u3[0][0] = u3_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u3[(r.size() - 1)][(z.size() - 1)] = u3_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the polar velocity is zero
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     rp = r[jr];
//     zp = z[0];
//     u3[jr][0] = u3_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u3[u1rSize - jr][0] = u3_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//     rp = r[jr];
//     zp = z[(z.size() - 1)];
//     u3[jr][(z.size() - 1)] = u3_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//     u3[u1rSize - jr][(z.size() - 1)] = u3_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       uz = sqrt(kerr_geo_Vz(a, En, Lz, Qc, zp))/(1. - zp*zp);
//       // take into account when the polar velocity is positive (z-velocity is negative)
//       u3[jr][jz] = u3_subfunction(a, En, Lz, Qc, rp, zp, -uz);
//       u3[u1rSize - jr][jz] = u3_subfunction(a, En, Lz, Qc, rp, zp, -uz);

//       // take into account when the polar velocity is negative (z-velocity is positive)
//       u3[jr][u1zSize - jz] = u3_subfunction(a, En, Lz, Qc, rp, zp, uz);
//       u3[u1rSize - jr][u1zSize - jz] = u3_subfunction(a, En, Lz, Qc, rp, zp, uz);
//     }
//   }

//   return u3;
// }

// ComplexMatrix tetrad_velocity_3_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u3(u1rSize, ComplexVector(u1zSize));
//   double rp, zp, uz, Vz, uzdz;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u3[0][0] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u3[(r.size() - 1)][(z.size() - 1)] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);

//   // take into account when the polar velocity is zero
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     rp = r[jr];
//     zp = z[0];
//     u3[jr][0] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);
//     u3[u1rSize - jr][0] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);

//     rp = r[jr];
//     zp = z[(z.size() - 1)];
//     u3[jr][(z.size() - 1)] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);
//     u3[u1rSize - jr][(z.size() - 1)] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];

//       Vz = kerr_geo_Vz(a, En, Lz, Qc, zp);
//       uz = sqrt(Vz)/(1. - zp*zp);
//       // Vzdz = kerr_geo_Vz_dz(a, En, Lz, Qc, zp);
//       // uzdz = (2.*zp/(1. - zp*zp) + 0.5*Vzdz/Vz)*uz;
//       uzdz = 0.; // based on choice of rigid extension of u^\alpha
//       // take into account when the polar velocity is positive (z-velocity is negative)
//       u3[jr][jz] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz);
//       u3[u1rSize - jr][jz] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz);

//       // take into account when the polar velocity is negative (z-velocity is positive)
//       u3[jr][u1zSize - jz] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz);
//       u3[u1rSize - jr][u1zSize - jz] = u3_dz_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz);
//     }
//   }

//   return u3;
// }

// ComplexMatrix tetrad_velocity_3_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u3(u1rSize, ComplexVector(u1zSize));
//   double rp, zp, uz, Vz, uzdz, uzdz2;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u3[0][0] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u3[(r.size() - 1)][(z.size() - 1)] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

//   // take into account when the polar velocity is zero
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     rp = r[jr];
//     zp = z[0];
//     u3[jr][0] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
//     u3[u1rSize - jr][0] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

//     rp = r[jr];
//     zp = z[(z.size() - 1)];
//     u3[jr][(z.size() - 1)] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
//     u3[u1rSize - jr][(z.size() - 1)] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];

//       Vz = kerr_geo_Vz(a, En, Lz, Qc, zp);
//       uz = sqrt(Vz)/(1. - zp*zp);
//       // Vzdz = kerr_geo_Vz_dz(a, En, Lz, Qc, zp);
//       // Vzdz2 = kerr_geo_Vz_dz2(a, En, Lz, Qc, zp);
//       // uzdz = (2.*zp/(1. - zp*zp) + 0.5*Vzdz/Vz)*uz;
//       // uzdz2 = (2.*(1 + 3.*zp*zp)/pow(1. - zp*zp, 2) + 2.*zp*Vzdz/Vz/(1. - zp*zp) - pow(0.5*Vzdz/Vz, 2) + 0.5*Vzdz2/Vz)*uz;
//       uzdz = 0.; // based on choice of rigid extension of u^\alpha
//       uzdz2 = 0.; // based on choice of rigid extension of u^\alpha
//       // take into account when the polar velocity is positive (z-velocity is negative)
//       u3[jr][jz] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz, -uzdz2);
//       u3[u1rSize - jr][jz] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz, -uzdz2);

//       // take into account when the polar velocity is negative (z-velocity is positive)
//       u3[jr][u1zSize - jz] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz, uzdz2);
//       u3[u1rSize - jr][u1zSize - jz] = u3_dz2_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz, uzdz2);
//     }
//   }

//   return u3;
// }

// ComplexMatrix tetrad_velocity_4(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u4(u1rSize, ComplexVector(u1zSize));
//   double rp, zp;
//   // first deal with double turning points
//   // rmin, zmin
//   rp = r[0];
//   zp = z[0];
//   u4[0][0] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // rmin, zmax
//   rp = r[0];
//   zp = z[(z.size() - 1)];
//   u4[0][(z.size() - 1)] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // rmax, zmax
//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u4[(r.size() - 1)][(z.size() - 1)] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // rmax, zmin
//   rp = r[(r.size() - 1)];
//   zp = z[0];
//   u4[(r.size() - 1)][0] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);

//   // take into account when the radial velocity is zero but different polar velocities
//   for(size_t jz = 1; jz < z.size() - 1; jz++){
//     rp = r[0];
//     zp = z[jz];
// 	double uz = sqrt(kerr_geo_Vz(a, En, Lz, Qc, zp))/(1. - zp*zp);
//     u4[0][jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, uz);
//     u4[0][u1zSize - jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, -uz);

//     rp = r[(r.size() - 1)];
//     zp = z[jz];
//     u4[(r.size() - 1)][jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, uz);
//     u4[(r.size() - 1)][u1zSize - jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, -uz);
//   }

//   // take into account when the polar velocity is zero but different radial velocities
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     rp = r[jr];
//     zp = z[0];
// 	// take into account when the radial velocity is positive
// 	u4[jr][0] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);
// 	u4[jr][z.size() - 1] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);

// 	// take into account when the radial velocity is negative
// 	u4[u1rSize - jr][0] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);
// 	u4[u1rSize - jr][z.size() - 1] = u4_subfunction(a, En, Lz, Qc, rp, zp, 0.);
//   }

//   // take into account different radial and polar velocities
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];
//       double uz = sqrt(kerr_geo_Vz(a, En, Lz, Qc, zp))/(1. - zp*zp);
//       // take into account when the radial velocity is positive
//       u4[jr][jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, uz);
//       u4[jr][u1zSize - jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, -uz);

//       // take into account when the radial velocity is negative
//       u4[u1rSize - jr][jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, uz);
//       u4[u1rSize - jr][u1zSize - jz] = u4_subfunction(a, En, Lz, Qc, rp, zp, -uz);
//     }
//   }

//   return u4;
// }

// ComplexMatrix tetrad_velocity_4_dz(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u4(u1rSize, ComplexVector(u1zSize));
//   double rp, zp, uz, Vz, uzdz;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u4[0][0] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u4[(r.size() - 1)][(z.size() - 1)] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);

//   // take into account when the polar velocity is zero
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     rp = r[jr];
//     zp = z[0];
//     u4[jr][0] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);
//     u4[u1rSize - jr][0] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);

//     rp = r[jr];
//     zp = z[(z.size() - 1)];
//     u4[jr][(z.size() - 1)] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);
//     u4[u1rSize - jr][(z.size() - 1)] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, 0., 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];

//       Vz = kerr_geo_Vz(a, En, Lz, Qc, zp);
//       uz = sqrt(Vz)/(1. - zp*zp);
//       // Vzdz = kerr_geo_Vz_dz(a, En, Lz, Qc, zp);
//       // uzdz = (2.*zp/(1. - zp*zp) + 0.5*Vzdz/Vz)*uz;
//       uzdz = 0.; // based on choice of rigid extension of u^\alpha
//       // take into account when the polar velocity is positive (z-velocity is negative)
//       u4[jr][jz] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz);
//       u4[u1rSize - jr][jz] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz);

//       // take into account when the polar velocity is negative (z-velocity is positive)
//       u4[jr][u1zSize - jz] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz);
//       u4[u1rSize - jr][u1zSize - jz] = u4_dz_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz);
//     }
//   }

//   return u4;
// }

// ComplexMatrix tetrad_velocity_4_dz2(double a, double En, double Lz, double Qc, Vector r, Vector z){
//   int u1rSize = 2*(r.size() - 1);
//   int u1zSize = 2*(z.size() - 1);
//   ComplexMatrix u4(u1rSize, ComplexVector(u1zSize));
//   double rp, zp, uz, Vz, uzdz, uzdz2;
//   // first deal with double turning points
//   rp = r[0];
//   zp = z[0];
//   u4[0][0] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

//   rp = r[(r.size() - 1)];
//   zp = z[(z.size() - 1)];
//   u4[(r.size() - 1)][(z.size() - 1)] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

//   // take into account when the polar velocity is zero
//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     rp = r[jr];
//     zp = z[0];
//     u4[jr][0] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
//     u4[u1rSize - jr][0] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);

//     rp = r[jr];
//     zp = z[(z.size() - 1)];
//     u4[jr][(z.size() - 1)] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
//     u4[u1rSize - jr][(z.size() - 1)] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, 0., 0., 0.);
//   }

//   for(size_t jr = 1; jr < r.size() - 1; jr++){
//     for(size_t jz = 1; jz < z.size() - 1; jz++){
//       rp = r[jr];
//       zp = z[jz];

//       Vz = kerr_geo_Vz(a, En, Lz, Qc, zp);
//       uz = sqrt(Vz)/(1. - zp*zp);
//       // Vzdz = kerr_geo_Vz_dz(a, En, Lz, Qc, zp);
//       // Vzdz2 = kerr_geo_Vz_dz2(a, En, Lz, Qc, zp);
//       // uzdz = (2.*zp/(1. - zp*zp) + 0.5*Vzdz/Vz)*uz;
//       // uzdz2 = (2.*(1 + 3.*zp*zp)/pow(1. - zp*zp, 2) + 2.*zp*Vzdz/Vz/(1. - zp*zp) - pow(0.5*Vzdz/Vz, 2) + 0.5*Vzdz2/Vz)*uz;
//       uzdz = 0.; // based on choice of rigid extension of u^\alpha
//       uzdz2 = 0.; // based on choice of rigid extension of u^\alpha
//       // take into account when the polar velocity is positive (z-velocity is negative)
//       u4[jr][jz] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz, -uzdz2);
//       u4[u1rSize - jr][jz] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, -uz, -uzdz, -uzdz2);

//       // take into account when the polar velocity is negative (z-velocity is positive)
//       u4[jr][u1zSize - jz] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz, uzdz2);
//       u4[u1rSize - jr][u1zSize - jz] = u4_dz2_subfunction(a, En, Lz, Qc, rp, zp, uz, uzdz, uzdz2);
//     }
//   }

//   return u4;
// }
