
//#include "fluxes.hpp"
//#include "waveform.hpp"
#include "inspiral.hpp"
#include "test.hpp"
#include "omp.h"
#include "ssf.hpp"
#include "unit_test.hpp"
#include "resflux.hpp"
#include "hertz.hpp"
#include "metric.hpp"
#include "redshift.hpp"

int main(){
	std::cout << std::setprecision(16);

	// GeodesicSource geo = kerr_geo_orbit(0.6, 8.284178275185655, 0.2, cos(0.25*M_PI), pow(2, 10));
	// int sampleN = pow(2, 6);
	// std::string dir = "../GEORG-Data/a6000_p8284_e2000_x0707_N64_LMAX35";
	// geo.save(dir);
	// save_regAB_data(geo, sampleN, dir);
	// save_params(35, sampleN, geo, dir);
	// ssf_components_parallel_m(geo, 35, sampleN, dir);
	// load_construct_and_save_ssf_multipoles(dir);

	// geo = kerr_geo_orbit(0.2, 10.041115133330957, 0.2, cos(0.25*M_PI), pow(2, 10));
	// sampleN = pow(2, 7);
	// dir = "../GEORG-Data/a2000_p10041_e2000_x0707_N256";
	// geo.save(dir);
	// save_regAB_data(geo, sampleN, dir);
	// save_params(25, sampleN, geo, dir);
	// ssf_components_parallel_m(geo, 25, sampleN, dir);
	// load_construct_and_save_ssf_multipoles(dir);
	//
	// GeodesicSource geo = kerr_geo_orbit(0.9, 4.607437338366917, 0.5, cos(0.25*M_PI), pow(2, 12));
	// int sampleN = pow(2, 7);
	// std::string dir = "../GEORG-Data/a9000_p4607_e5000_x7071_N128_v2";
	// geo.save(dir);
	// save_regAB_data(geo, sampleN, dir);
	// save_params(25, sampleN, geo, dir);
	// ssf_components_parallel_m(geo, 25, sampleN, dir);
	// load_construct_and_save_ssf_multipoles(dir);
	//
	// geo = kerr_geo_orbit(0.9, 7.787543300602265, 0.2, cos(M_PI/3.), pow(2, 10));
	// sampleN = pow(2, 7);
	// dir = "../GEORG-Data/a9000_p7788_e2000_x0500_N256";
	// geo.save(dir);
	// save_regAB_data(geo, sampleN, dir);
	// save_params(25, sampleN, geo, dir);
	// ssf_components_parallel_m(geo, 25, sampleN, dir);
	// load_construct_and_save_ssf_multipoles(dir);
	//
	// generate_adiabatic_inspiral_data_2d("../GEORG-Data/inspiral_newParameters_v5/", 65, 513);

	GeodesicSource geo = kerr_geo_orbit(0., 20., 0.3, 1, pow(2, 10));
	// ComplexVector res = resonant_flux_mode_coefficients(-2, 10, 10, 2, 3, 2, 3, geo);
	// Fluxes eflux2 = energy_flux_lm(-2, 2, 2, geo);
	// Fluxes efluxm2 = energy_flux_lm(-2, 2, -2, geo);
	// std::cout << eflux2.infinity + efluxm2.infinity << "\n";
	// ResonantFluxList fluxes = res_flux_lm(0, 2, 2, 3, 2, geo);
	// for(int i = 0; i < 9; i++){
	// 	std::cout << fluxes.Qdot[i].infinity << "\n";
	// }
	// res_flux_parallel_l(0, 3, 2, geo, 25, "../GEORG-Data/flux_test/");
	// full_flux_parallel_l(-2, geo, 30, "../GEORG-Data/flux_test/");
	// FluxList fluxes2 = flux_lm(0, 8, 8, geo);
	// std::cout << fluxes2.Edot.infinity << "\n";

	auto t1 = std::chrono::high_resolution_clock::now();
	RadialTeukolsky Rt(0., -2, 60, 2, 0.5, geo.getRadialPosition());
	Rt.generateSolutions(AUTO);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "Teukolsky took "
	          << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	          << " milliseconds\n";

	// int s = -2;
	// int L = 10;
	// int m = 10;
	// int k = 2;
	// int n = 3;
	// RadialTeukolsky Rt(geo.getBlackHoleSpin(), s, L, m, geo.getTimeFrequency(m, k, n), geo.getRadialPosition());
	// Rt.generateSolutions(AUTO);
	// SpinWeightedHarmonic Slm(s, L, m, geo.getBlackHoleSpin()*geo.getTimeFrequency(m, k, n), geo.getPolarPosition());
	// Slm.generateSolutions();
	// TeukolskyMode teuk(s, L, m, k, n, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// teuk.generateSolutions(Slm, Rt, geo);
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";

	// GeodesicSource geo;
	//
	// GeodesicSource geoCirc = kerr_geo_orbit(0.9, 10., 0., 1., pow(2, 2));
	//
	// RealMatrix xp = worldline_grid(0., 9., 11., 202, -0.99, 0.99, 21, 0., geoCirc);
	// Vector t = xp[0];
	// Vector r = xp[1];
	// Vector z = xp[2];
	// Vector phi = xp[3];
	// RealMatrix hmunu = metric_perturbation_circ(ORG, 2, 2, t, r, z, phi, geoCirc);
	// save_full_metric_data("data/hab_kerr_ORG.txt", t, r, z, phi, hmunu);
	// hmunu = metric_perturbation_circ(IRG, 2, 2, t, r, z, phi, geoCirc);
	// save_full_metric_data("data/hab_kerr_IRG.txt", t, r, z, phi, hmunu);
	// hmunu = metric_perturbation_circ(SAAB0, 2, 2, t, r, z, phi, geoCirc);
	// save_full_metric_data("data/hab_kerr_SAAB0.txt", t, r, z, phi, hmunu);
	// hmunu = metric_perturbation_circ(SAAB4, 2, 2, t, r, z, phi, geoCirc);
	// save_full_metric_data("data/hab_kerr_SAAB4.txt", t, r, z, phi, hmunu);
	// hmunu = metric_perturbation_circ(ASAAB0, 2, 2, t, r, z, phi, geoCirc);
	// save_full_metric_data("data/hab_kerr_ASAAB0.txt", t, r, z, phi, hmunu);
	// hmunu = metric_perturbation_circ(ASAAB4, 2, 2, t, r, z, phi, geoCirc);
	// save_full_metric_data("data/hab_kerr_ASAAB4.txt", t, r, z, phi, hmunu);

	// int lmax = 15;
	// std::string filepath = "data/huu_kerr/h33_test_kerr";
	// redshift_circular(filepath, SAAB, lmax, geoCirc);
	// redshift_circular(filepath, ASAAB, lmax, geoCirc);
	// redshift_circular(filepath, ORG, lmax, geoCirc);
	// redshift_circular(filepath, IRG, lmax, geoCirc);
	// std::cout << geoCirc.getTimeFrequency(3) << "\n";

	// Complex ampI, ampH;

	// std::cout << "(0, 10, 0, 1)\n";
	// ampI = -0.096869111783302688482248 + 0.034691088866233852517826*I;
	// ampH = 0.00085249088747438226386116 - 0.00021231389467651976331909*I;
	// GeodesicSource geo = kerr_geo_orbit(0., 10., 0., 1., pow(2, 4));
	// TeukolskyMode teuk(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// // std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// // std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// // std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// // std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0, 10.0, 0, Cos[\[Pi]/4.])\n";
	// ampI = -0.070574319833484021478549 + 0.025274310416867242591122*I;
	// ampH = 0.00062108512651933021254542 - 0.00015468200783664169003028*I;
	// geo = kerr_geo_orbit(1.e-14, 10., 0., cos(M_PI/4.), pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0, 10.0, 0.1, 1)\n";
	// ampH = 0.00078064421192569808851129185 - 0.00019244803540255800765595629*I;
	// ampI = -0.090966715053256986219528369 + 0.032357265869968441307292536*I;
	// geo = kerr_geo_orbit(0., 10., 0.1, 1., pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0, 10.0, 0.1, Cos[\[Pi]/4.])\n";
	// ampI = -0.066274108683182594817608793 + 0.023573995759893462190169796*I;
	// ampH = 0.00056874098744548149271895618 - 0.00014020866870554800968163363*I;
	// geo = kerr_geo_orbit(1.e-14, 10., 0.1, cos(M_PI/4.), pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0.1, 10.0, 0, 1)\n";
	// ampI = -0.096227584894766381938308 + 0.034422505898013111362993*I;
	// ampH = 0.00082969529874874947045813 - 0.00004165341379611482771148*I;
	// geo = kerr_geo_orbit(0.1, 10., 0., 1., pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0.1, 10.0, 0, Cos[\[Pi]/4.]])\n";
	// ampI = -0.070468528422932792881 + 0.025252522903963179701*I;
	// ampH = 0.00060463777485046264504 - 0.00003075360976614048730*I;
	// geo = kerr_geo_orbit(0.1, 10., 0., cos(M_PI/4.), pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0.1, 10.0, 0, Cos[\[Pi]/4.]])\n";
	// ampI = 2.817674016780767525e-6 - 1.506523056130186854e-6*I;
	// ampH = -9.033598925366890723e-9 + 2.194950369301104043e-9*I;
	// geo = kerr_geo_orbit(0.1, 10., 0., cos(M_PI/4.), pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 2, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0.1, 10.0, 0.1, 1)\n";
	// ampH = 0.00076250722962311750779893683 - 0.00003624478217354641654074258*I;
	// ampI = -0.090689967567633577108184261 + 0.032215375076956083066057133*I;
	// geo = kerr_geo_orbit(0.1, 10., 0.1, 1., pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << "(0.1, 10.0, 0.1, Cos[\[Pi]/4.]])\n";
	// ampI = -0.066318324275833298296962 + 0.0236010362126117048938706*I;
	// ampH = 0.000554901107100978524115932 - 0.000026750779539412935880262*I;
	// geo = kerr_geo_orbit(0.1, 10., 0.1, cos(M_PI/4.), pow(2, 8));
	// teuk = TeukolskyMode(0, 2, 2, 0, 0, geo);
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(In)/ampH) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << abs(1. - teuk.getTeukolskyAmplitude(Up)/ampI) << "\n";
	//
	// std::cout << teuk_in_static_solution(0.9, 2, 3, 0, 10.) << "\n";
	// std::cout << teuk_in_static_derivative(0.9, 2, 3, 0, 10.) << "\n";
	// std::cout << teuk_up_static_solution(0.9, 2, 3, 0, 10.) << "\n";
	// std::cout << teuk_up_static_derivative(0.9, 2, 3, 0, 10.) << "\n";
	// std::cout << teuk_in_static_solution(0.9, -2, 3, 0, 10.) << "\n";
	// std::cout << teuk_in_static_derivative(0.9, -2, 3, 0, 10.) << "\n";
	// std::cout << teuk_up_static_solution(0.9, -2, 3, 0, 10.) << "\n";
	// std::cout << teuk_up_static_derivative(0.9, -2, 3, 0, 10.) << "\n";


	// GeodesicSource geo = kerr_geo_orbit(0.9, 8.1, 0.4, 0.5, pow(2, 11));
	// TeukolskyMode teuk(-2, 2, 1, 0, 0, geo);
	// // int sample = 160;
	// teuk.generateSolutions(geo);
	// std::cout << teuk.getTeukolskyAmplitude(In) << "\n";
	// std::cout << teuk.getTeukolskyAmplitude(Up) << "\n";
	// std::cout << "Teukolsky mode generated\n";
	// HertzMode hertz(teuk);
	// hertz.generateSolutions();
	// std::cout << hertz.getHertzAmplitude(In) << "\n";
	// std::cout << hertz.getHertzAmplitude(Up) << "\n";
	// std::cout << "Hertz mode generated\n";
	// ComplexTensor huuCoeff = redshift_coefficients(geo, pow(2, 7));
	// std::cout << "Redshift coefficients generated\n";
	//
	// SphericalHarmonicCoupling Cjlm(25, 2);
	// Cjlm.generateCouplings();
	//
	// auto t1 = std::chrono::high_resolution_clock::now();
	// for(size_t jr = 0; jr < huuCoeff[0].size()/2 + 1; jr++){
	// 	for(size_t jz = 0; jz < huuCoeff[0].size()/2 + 1; jz++){
	// 		for(int sgnur = 1; sgnur >= -1; sgnur-=2){
	// 			for(int sgnuz = 1; sgnuz >= -1; sgnuz-=2){
	// 				// std::cout << "huu(jr = "<<jr<<", jz = "<<jz<<"; sgnur = "<<sgnur<<", sgnuz = "<<sgnuz<<") = " << redshift_mode(2, hertz.getAzimuthalModeNumber(), hertz.getFrequency(), jr, jz, sgnur, sgnuz, huuCoeff, hertz, In) << "\n";
	// 				for(int l = 2; l <= 25; l++){
	// 					redshift_mode(huuCoeff, Cjlm, hertz, In, 2, jr, jz, sgnur, sgnuz);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
  // auto t2 = std::chrono::high_resolution_clock::now();
  // std::cout << "f() took "
  //           << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
  //           << " milliseconds\n";

	return 0;
}
