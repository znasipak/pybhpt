#ifndef SSF_HPP
#define SSF_HPP

#include "boost/filesystem.hpp"
#include "teukolsky.hpp"
#include "regularization.hpp"

typedef struct SelfForceDataStruct{
  RealTensor in;
  RealTensor up;
} SelfForceData;

typedef struct ParamsFileContentsStruct{
	int lmax;
	int orbitSampleNum;
	int sfSampleNum;
	double a;
	double p;
	double e;
	double x;
} ParamsFileContents;

void polar_coupling_coefficient(int order, Vector &couplingVector, TeukolskyMode &teuk, int l, int m);
double polar_coupling_coefficient(int order, TeukolskyMode &teuk, int l, int m, int jth);

double beta_coupling_higher_order(int n, int l, int m, double theta);
double beta_coupling(int n, int l, int m, double theta);

Vector scalar_self_force_reference(int mu, GeodesicSource &geo,
  int sampleSize);
RealTensor scalar_self_force_reference(List components, int lmax, int m,
  GeodesicSource &geo, int sampleSize);

RealTensor scalar_self_force_components_data_init(int componentNum, int lmax,
  int m, int sampleSize);
RealTensor scalar_self_force_components_convergence_init(const RealTensor &convergenceData, int convergenceCriteria);

int scalar_self_force_components_convergence_sum(RealTensor &ssfTot,
  RealTensor &ssfMode, RealTensor &ssfPrevious, RealTensor &ssfRef, RealTensor &ssfConvergence, int convergenceCriteria);
int scalar_self_force_components_convergence_sum(RealMatrix &ssfTot,
  RealMatrix &ssfMode, RealMatrix &ssfPrevious, RealMatrix &ssfRef, RealMatrix &ssfConvergence, int convergenceCriteria);

SelfForceData scalar_self_force_components_m(int mu, int lmax, int m,
  GeodesicSource &geo, int sampleSize);
SelfForceData scalar_self_force_components_m(List components, int lmax, int m,
  GeodesicSource &geo, int sampleSize);
int scalar_self_force_components_m(List components, RealTensor &ssfTotIn,
  RealTensor &ssfTotUp, RealTensor &ssfRef, int m, GeodesicSource geo);

SelfForceData scalar_self_force_components_mk(List components, int lmax, int m, int k, GeodesicSource &geo, int sampleSize);
int scalar_self_force_components_mk(List components, RealTensor &ssfTotIn, RealTensor &ssfTotUp, const RealTensor &ssfConvergenceInM, const RealTensor &ssfConvergenceUpM, int convergenceCriteriaM, RealTensor &ssfRef, int m, int k, GeodesicSource &geo);

double beta_coupling(int coeff, int l, int m, double thp);
double polar_coupling_coefficient(TeukolskyMode &teuk, int l, int m, int jth);
Complex self_force_amplitude(int mu, TeukolskyMode &teuk, int l, int m, int jr,
  int jth, BoundaryCondition bc);
RealMatrix scalar_self_force_components_mkn(int mu, int lmax, int m, int k,
  int n, GeodesicSource &geo, int sampleSize);
SelfForceData scalar_self_force_components_mkn(List components, int lmax,
  int m, int k, int n, GeodesicSource &geo, int sampleSize);


int scalar_self_force_components_mkn(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_equatorial(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_spherical(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_circular(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_generic(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, int m, int k, int n, GeodesicSource &geo);

int scalar_self_force_components_mkn(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp,
  int convergenceCriteria, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_equatorial(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp,
  int convergenceCriteria, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_spherical(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp,
  int convergenceCriteria, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_circular(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp,
  int convergenceCriteria, int m, int k, int n, GeodesicSource &geo);
int scalar_self_force_components_mkn_generic(List components, RealTensor &ssfIn,
  RealTensor &ssfUp, const RealTensor &ssfConvergenceIn, const RealTensor &ssfConvergenceUp,
  int convergenceCriteria, int m, int k, int n, GeodesicSource &geo);

int save_ssf_data_lm_ecceq(ComplexVector ssfData, int l, int m, GeodesicSource geo, const std::string &dir);

int save_ssf_data(List components, SelfForceData ssfData, GeodesicSource geo, const std::string &dir);
int save_ssf_data(List components, SelfForceData ssfData, int m, GeodesicSource geo, const std::string &dir);
int save_ssf_data(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag = 1);
int save_ssf_data_equatorial(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag = 1);
int save_ssf_data_spherical(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag = 1);
int save_ssf_data_circular(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag = 1);
int save_ssf_data_generic(List components, RealTensor ssfData, int m, GeodesicSource geo, const std::string &dir, int lmFlag = 1);

int save_params(int lmax, int sampleNum, GeodesicSource& geo, const std::string &dir);
ParamsFileContents load_params_file(const std::string &dir);

SelfForceData load_ssf_data(const std::string &dir);
SelfForceData load_ssf_data(int m, const std::string &dir);
SelfForceData load_ssf_data(int l, int m, const std::string &dir);

SelfForceData load_and_construct_ssf_multipoles(const std::string &dir);
int load_construct_and_save_ssf_multipoles(const std::string &dir);

#endif
