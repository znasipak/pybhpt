//unit_test.cpp

#include "unit_test.hpp"

void run_unit_tests(){
  int start, stop;
  start = clock();
  std::cout << "(UNIT) Running series of unit tests...\n";
  int testSuccessCount = 0;
  int testNum = 0;
  if(test_kerr_circ_eq_flux() > 0){
    testSuccessCount += 1;
    testNum += 1;
  }else{
    testNum += 1;
  }
  if(test_kerr_circ_eq_flux_2() > 0){
    testSuccessCount += 1;
    testNum += 1;
  }else{
    testNum += 1;
  }
  if(test_schwarzschild_circ_eq_flux() > 0){
    testSuccessCount += 1;
    testNum += 1;
  }else{
    testNum += 1;
  }
  if(test_teukolsky_radial_function() > 0){
    testSuccessCount += 1;
    testNum += 1;
  }else{
    testNum += 1;
  }
  if(test_teukolsky_radial_function_highModes() > 0){
    testSuccessCount += 1;
    testNum += 1;
  }else{
    testNum += 1;
  }
  if(test_teukolsky_polar_function() > 0){
    testSuccessCount += 1;
    testNum += 1;
  }else{
    testNum += 1;
  }
  stop = clock();
  double duration = (stop-start)/double(CLOCKS_PER_SEC);
  std::cout << "(UNIT) " << testSuccessCount << " out of " << testNum << " unit tests passed.\n";
  std::cout << "(UNIT) Tests completed in "<<duration<<" seconds. \n";
}

int test_kerr_circ_eq_flux(){
  int test = 0;
  double tol = 1.e-9;

  double testOrbitSpin = 0.9;
  double testOrbitRadius = 3.237639327804;
  double testValue_EdotInfinity = 0.00907246856435667;
  double testValue_EdotHorizon = -0.0002053208822379471;

  GeodesicSource geo = kerr_geo_circ(testOrbitSpin, testOrbitRadius, 1);
  Fluxes edot = energy_flux(-2, geo);
  if(abs(1. - testValue_EdotInfinity/edot.infinity) < tol){
    test += 1;
  }
  if(abs(1. - testValue_EdotHorizon/edot.horizon) < tol){
    test += 1;
  }

  int testTot = 2;
  if(test == testTot){
    std::cout << "(UNIT) PASSED: Energy flux for (a, p, e, x) = ("<<testOrbitSpin<<", "<<testOrbitRadius<<", 0, 1). \n";
  }else{
    std::cout << "(UNIT) FAILED: Energy flux for (a, p, e, x) = ("<<testOrbitSpin<<", "<<testOrbitRadius<<", 0, 1). \n";
  }

  return test;
}

int test_kerr_circ_eq_flux_2(){
  int test = 0;
  double tol = 1.e-9;

  double testOrbitSpin = 0.9;
  double testOrbitRadius = 6.687053095016;
  double testValue_EdotInfinity = 0.00034094275249432537;
  double testValue_EdotHorizon = -2.0087269373724617e-6;

  GeodesicSource geo = kerr_geo_circ(testOrbitSpin, testOrbitRadius, 1);
  Fluxes edot = energy_flux(-2, geo);
  if(abs(1. - testValue_EdotInfinity/edot.infinity) < tol){
    test += 1;
  }
  if(abs(1. - testValue_EdotHorizon/edot.horizon) < tol){
    test += 1;
  }

  int testTot = 2;
  if(test == testTot){
    std::cout << "(UNIT) PASSED: Energy flux for (a, p, e, x) = ("<<testOrbitSpin<<", "<<testOrbitRadius<<", 0, 1). \n";
  }else{
    std::cout << "(UNIT) FAILED: Energy flux for (a, p, e, x) = ("<<testOrbitSpin<<", "<<testOrbitRadius<<", 0, 1). \n";
  }

  return test;
}

int test_schwarzschild_circ_eq_flux(){
  int test = 0;
  double tol = 1.e-9;
  double testOrbitSpin = 0.;
  double testOrbitRadius = 6.105567305770028;
  double testValue_EdotInfinity = 0.0008492678896009272;
  double testValue_EdotHorizon = 2.5004530589418276e-6;

  GeodesicSource geo = kerr_geo_circ(testOrbitSpin, testOrbitRadius, 1);
  Fluxes edot = energy_flux(-2, geo);
  if(abs(1. - testValue_EdotInfinity/edot.infinity) < tol){
    test += 1;
  }
  if(abs(1. - testValue_EdotHorizon/edot.horizon) < tol){
    test += 1;
  }

  int testTot = 2;
  if(test == testTot){
    std::cout << "(UNIT) PASSED: Energy flux for (a, p, e, x) = ("<<testOrbitSpin<<", "<<testOrbitRadius<<", 0, 1). \n";
  }else{
    std::cout << "(UNIT) FAILED: Energy flux for (a, p, e, x) = ("<<testOrbitSpin<<", "<<testOrbitRadius<<", 0, 1). \n";
  }

  return test;
}

int test_teukolsky_radial_function(){
  int test = 0;
  double tol = 1.e-11;
  int s = -2;
  int l = 2;
  int m = 2;
  double a = 0.9;
  double om = 0.15;

  Vector r = {4., 10.};
  RadialTeukolsky teuk(a, s, l, m, om, r);
  teuk.generateSolutions();
  Complex RinMin = teuk.getSolution(In, 0);
  Complex RinMax = teuk.getSolution(In, 1);
  Complex RupMin = teuk.getSolution(Up, 0);
  Complex RupMax = teuk.getSolution(Up, 1);
  Complex RinMinP = teuk.getDerivative(In, 0);
  Complex RinMaxP = teuk.getDerivative(In, 1);
  Complex RupMinP = teuk.getDerivative(Up, 0);
  Complex RupMaxP = teuk.getDerivative(Up, 1);

  Complex testValue_RinMin = 83.15413293250249 - 22.88043485716941*I;
  Complex testValue_RinMax = 4755.347134819352 + 4033.570104491953*I;
  Complex testValue_RupMin = 534.2454295305542 - 27.8745682775869*I;
  Complex testValue_RupMax = -439.0287646624466 - 144.4480522305074*I;
  Complex testValue_RinMinP = 118.36072943717715 - 8.51321106827702*I;
  Complex testValue_RinMaxP = 1448.239361624757 + 2280.937411825064*I;
  Complex testValue_RupMinP = -182.8647917286560 - 129.1459266988753*I;
  Complex testValue_RupMaxP = -252.8589469480107 - 37.6026300231709*I;

  if(abs(1. - RinMin/testValue_RinMin) < tol){
    test += 1;
  }
  if(abs(1. - RinMax/testValue_RinMax) < tol){
    test += 1;
  }
  if(abs(1. - RupMin/testValue_RupMin) < tol){
    test += 1;
  }
  if(abs(1. - RupMax/testValue_RupMax) < tol){
    test += 1;
  }
  if(abs(1. - RinMinP/testValue_RinMinP) < tol){
    test += 1;
  }
  if(abs(1. - RinMaxP/testValue_RinMaxP) < tol){
    test += 1;
  }
  if(abs(1. - RupMinP/testValue_RupMinP) < tol){
    test += 1;
  }
  if(abs(1. - RupMaxP/testValue_RupMaxP) < tol){
    test += 1;
  }

  int testTot = 8;
  if(test == testTot){
    std::cout << "(UNIT) PASSED: Radial Teukolsky solutions for parameters (a, s, l, m, omega) = ("<<a<<", "<<s<<", "<<l<<", "<<m<<", "<<om<<")\n";
  }else{
    std::cout << "(UNIT) FAILED: Radial Teukolsky solutions for parameters (a, s, l, m, omega) = ("<<a<<", "<<s<<", "<<l<<", "<<m<<", "<<om<<"). Only passed " << test << " of " << testTot << " tests \n";
  }

  return test;
}

int test_teukolsky_radial_function_highModes(){
  int test = 0;
  double tol = 1.e-11;
  int s = -2;
  int l = 12;
  int m = 2;
  double a = 0.9;
  double om = 0.65;

  Vector r = {4., 10.};
  RadialTeukolsky teuk(a, s, l, m, om, r);
  teuk.generateSolutions();
  Complex RinMin = teuk.getSolution(In, 0);
  Complex RinMax = teuk.getSolution(In, 1);
  Complex RupMin = teuk.getSolution(Up, 0);
  Complex RupMax = teuk.getSolution(Up, 1);
  Complex RinMinP = teuk.getDerivative(In, 0);
  Complex RinMaxP = teuk.getDerivative(In, 1);
  Complex RupMinP = teuk.getDerivative(Up, 0);
  Complex RupMaxP = teuk.getDerivative(Up, 1);

  Complex testValue_RinMin = Complex(2.405458988239966e11,2.080764223659599e11);
  Complex testValue_RinMax = Complex(9.69925822204652e16,5.276346454776849e17);
  Complex testValue_RupMin = Complex(1.282287732253362e9,-3.245083146650749e9);
  Complex testValue_RupMax = Complex(-25190.65925803417,-62249.83277489972);
  Complex testValue_RinMinP = Complex(1.076190957914128e12,9.78353779646494e11);
  Complex testValue_RinMaxP = Complex(6.58615110099728e16,7.109512466221851e17);
  Complex testValue_RupMinP = Complex(-4.929425935275619e9,1.1315542345652564e10);
  Complex testValue_RupMaxP = Complex(15136.92639217653,62776.41393296017);

  // std::cout << abs(1. - RinMin/testValue_RinMin) << "\n";
  // std::cout << abs(1. - RinMax/testValue_RinMax) << "\n";
  // std::cout << abs(1. - RupMin/testValue_RupMin) << "\n";
  // std::cout << abs(1. - RupMax/testValue_RupMax) << "\n";

  if(abs(1. - RinMin/testValue_RinMin) < tol){
    test += 1;
  }
  if(abs(1. - RinMax/testValue_RinMax) < tol){
    test += 1;
  }
  if(abs(1. - RupMin/testValue_RupMin) < tol){
    test += 1;
  }
  if(abs(1. - RupMax/testValue_RupMax) < tol){
    test += 1;
  }
  if(abs(1. - RinMinP/testValue_RinMinP) < tol){
    test += 1;
  }
  if(abs(1. - RinMaxP/testValue_RinMaxP) < tol){
    test += 1;
  }
  if(abs(1. - RupMinP/testValue_RupMinP) < tol){
    test += 1;
  }
  if(abs(1. - RupMaxP/testValue_RupMaxP) < tol){
    test += 1;
  }

  int testTot = 8;
  if(test == testTot){
    std::cout << "(UNIT) PASSED: Radial Teukolsky solutions for parameters (a, s, l, m, omega) = ("<<a<<", "<<s<<", "<<l<<", "<<m<<", "<<om<<") \n";
  }else{
    std::cout << "(UNIT) FAILED: Radial Teukolsky solutions for parameters (a, s, l, m, omega) = ("<<a<<", "<<s<<", "<<l<<", "<<m<<", "<<om<<"). Only passed " << test << " of " << testTot << " tests \n";
  }

  return test;
}

int test_teukolsky_polar_function(){
  int test = 0;
  double tol = 1.e-11;
  int s = 1;
  int l = 3;
  int m = 2;
  double a = 0.9;
  double om = 0.65;

  Vector th = {.1, 1.};
  SpinWeightedHarmonic Slm(s, l, m, a*om, th);
  Slm.generateSolutionsAndDerivatives();
  double SlmMin = Slm.getSolution(0);
  double SlmMax = Slm.getSolution(1);
  double SlmMinP = Slm.getDerivative(0);
  double SlmMaxP = Slm.getDerivative(1);

  double testValue_SlmMin = 0.0005140627326102776;
  double testValue_SlmMax = 0.2824947558777588;
  double testValue_SlmMinP = 0.015364831841802413;
  double testValue_SlmMaxP = 0.47201677488912486;

  if(abs(1. - SlmMin/testValue_SlmMin) < tol){
    test += 1;
  }
  if(abs(1. - SlmMax/testValue_SlmMax) < tol){
    test += 1;
  }
  if(abs(1. - SlmMinP/testValue_SlmMinP) < tol){
    test += 1;
  }
  if(abs(1. - SlmMaxP/testValue_SlmMaxP) < tol){
    test += 1;
  }

  int testTot = 4;
  if(test == testTot){
    std::cout << "(UNIT) PASSED: Radial Teukolsky solutions for parameters (a, s, l, m, omega) = ("<<a<<", "<<s<<", "<<l<<", "<<m<<", "<<om<<") \n";
  }else{
    std::cout << "(UNIT) FAILED: Radial Teukolsky solutions for parameters (a, s, l, m, omega) = ("<<a<<", "<<s<<", "<<l<<", "<<m<<", "<<om<<"). Only passed " << test << " of " << testTot << " tests \n";
  }

  return test;
}
