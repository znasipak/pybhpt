//#include "fluxes.hpp"
//#include "waveform.hpp"
#include "test.hpp"

void test_special_functions(){
	Complex a = 50.2 + 5.76*I;
	Complex b = -24.2 - 12.*I;
	Complex z = -3.45 + 4.1*I;
	std::cout << std::setprecision(15);

	Result U = hyper_u(a, b, z);
	std::cout << "U(a, b, z) = " << U << "\n";
	Result U2 = hyper_u(a, b, -z);
	std::cout << "U(a, b, -z) = " << U2 << "\n";
	Result Ai = hyper_Ai(z);
	std::cout << "Ai(z) = " << Ai << "\n";
	Result J = hyper_bessel_J(a, z);
	std::cout << "J_a(z) = " << J << "\n";
	Result K = hyper_bessel_K(b, z);
	std::cout << "K_b(z) = " << K << "\n";
	Result G = hyper_incomplete_Gamma(b, z);
	std::cout << "G_b(z) = " << G << "\n";
	Result erfN = hyper_erf(z);
	std::cout << "erf(z) = " << erfN << "\n";
	Result erfcN = hyper_erfc(z);
	std::cout << "erfc(z) = " << erfcN << "\n";
}

void ssf_components_parallel_m(GeodesicSource geo, int lmax, int mminTemp, int mmaxTemp, int sampleSSF, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
		boost::filesystem::create_directory(dir);
	}

	// int modeNumMax = (lMax - 1)*(lMax + 4)/2;
	// int modeMax = 2*lmax + 1;
	// std::vector<SelfForceData> ssfMmodes(modeMax);
	int m;

	std::string dirUp;
	std::string dirIn;
	std::string dirReg;
	if(dir.back() == '/'){
		dirUp = dir + "up";
		dirIn = dir + "in";
		dirReg = dir + "reg";
	}else{
		dirUp = dir + "/up";
		dirIn = dir + "/in";
		dirReg = dir + "/reg";
	}

	int mmax = mmaxTemp;
	int mmin = mminTemp;
	if(mmax < mmin){
		mmax = mmin;
		mmin = mmaxTemp;
	}
	if(mmax > lmax){
		mmax = lmax;
	}
	if(mmin < -lmax){
		mmin = -lmax;
	}

	save_regAB_data(geo, sampleSSF, dirReg);
	// #pragma omp parallel shared(geo, ssfMmodes, modeMax, lmax) private(m)
	#pragma omp parallel num_threads(15) shared(geo, mmax, mmin, sampleSSF) private(m)
	{
		#pragma omp for schedule(dynamic)
			for(m = mmax; m >= mmin; m--) {
				int th_id = omp_get_thread_num();
				sleep(abs(m));
				std::cout << "Start SSF calculation for ("<<m<<")-mode on kernel " << th_id << "\n";
				// ssfMmodes[lmax + m] = scalar_self_force_components_m(components, lmax, m, geo, sampleSSF);
				// std::cout << "SSF calculation for ("<<m<<")-mode completed on kernel " << th_id << "\n";
				// int status = 0;
				// status = save_ssf_data(components, ssfMmodes[lmax + m].up, m, geo, dirUp);
				// status = save_ssf_data(components, ssfMmodes[lmax + m].in, m, geo, dirIn);
				List components = {0, 1, 2, 3, 4};
				// List components = {3, 0, 1, 2};
				// List components = {0, 3};
				SelfForceData ssfMmode = scalar_self_force_components_m(components, lmax, m, geo, sampleSSF);
				std::cout << "SSF calculation for ("<<m<<")-mode completed on kernel " << th_id << "\n";
				int status = 0;
				status = save_ssf_data(components, ssfMmode.up, m, geo, dirUp);
				status = save_ssf_data(components, ssfMmode.in, m, geo, dirIn);
			}
	}
}

void ssf_components_parallel_m(GeodesicSource geo, int lmax, int sampleSSF, const std::string& dir){
	ssf_components_parallel_m(geo, lmax, -lmax, lmax, sampleSSF, dir);
}

void ssf_components_m(GeodesicSource geo, int lmax, int m, int sampleSSF, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
		boost::filesystem::create_directory(dir);
	}

	List components = {0, 1, 2, 3, 4};
	// List components = {1};
	std::string dirUp = dir + "/up";
	std::string dirIn = dir + "/in";

	std::cout << "Start SSF calculation for ("<<m<<")-mode on main kernel \n";
	SelfForceData ssf = scalar_self_force_components_m(components, lmax, m, geo, sampleSSF);
	std::cout << "SSF calculation for ("<<m<<")-mode completed on main kernel \n";
	int status = 0;
	status = save_ssf_data(components, ssf.up, m, geo, dirUp);
	status = save_ssf_data(components, ssf.in, m, geo, dirIn);
}

void test_teukolsky_ZlmUp_PN(const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	int lmin = 2;
	int lmax = 3;
	int mmin = 2;
	int nmin = 1;
	int nmax = 2;

	// double lemin = log(1e-3);
	// double lemax = log(0.1);
	// double ecount = 4;
	// double leinc = (lemax - lemin)/(ecount - 1.);
	//
	// double lpmin = log(200);
	// double lpmax = log(2e4);
	// double pcount = 10;
	// double lpinc = (lpmax - lpmin)/(pcount - 1.);

	double a = 0.9;
	double x = 1.;

	double eList[4] = {0.01, 0.1, 0.25, 0.5};
	double pList[2] = {10, 20};

	int sampleSize = pow(2, 9);

	// RealTensor ReZlmn(lmax - lmin + 1, RealMatrix(lmax + 1, Vector(nmax - nmin + 1, 0.)));
	// RealTensor ImZlmn = ReZlmn;
	Complex Zlmn = 0.;

	for(int l = lmin; l <= lmax; l++){
		for(int m = mmin; m <= l; m++){
			for(int n = nmin; n <= nmax; n++){
				char buff[9];
				std::snprintf(buff, 9, "%d_%d_%d", l, m, n);
				std::string filename = "ZlmnUp_";
				filename += buff;
				filename += ".txt";
				std::string filepath = dir + "/" + filename;

				std::ofstream file;
				file.open(filepath);

				file << "a\tl\tm\tn\n";
				file << a << "\t" << l << "\t" << m << "\t" << n << "\n\n";
				file << "p\te\tRe[Zlmn^+]\tIm[Zlmn^+]\n";

				// for(double lp = lpmin; lp <= lpmax; lp += lpinc){
				// 	for(double le = lemin; le <= lemax; le += leinc){
				// 		std::cout << "Generating orbit for (a, p, e, x) = ("<<a<<", "<<exp(lp)<<", "<<exp(le)<<", "<<x<<") \n";
				// 		GeodesicSource geo = kerr_geo_orbit(a, exp(lp), exp(le), x, sampleSize);
				// 		std::cout << "Generating Teukolsky amplitude for (l, m, n) = ("<<l<<", "<<m<<", "<<n<<") \n";
				// 		TeukolskyMode teuk(-2, l, m, 0, n, geo);
				// 		teuk.generateSolutions(geo);
				// 		Zlmn = teuk.getTeukolskyAmplitude(Up);
				// 		file << std::scientific << std::setprecision(15);
				// 		file << exp(lp) << "\t" << exp(le) << "\t" << std::real(Zlmn) <<"\t" << std::imag(Zlmn) << "\n";
				// 	}
				// }

				for(double p : pList){
					for(double e : eList){
						std::cout << "Generating orbit for (a, p, e, x) = ("<<a<<", "<<p<<", "<<e<<", "<<x<<") \n";
						GeodesicSource geo = kerr_geo_orbit(a, p, e, x, sampleSize);
						std::cout << "Generating Teukolsky amplitude for (l, m, n) = ("<<l<<", "<<m<<", "<<n<<") \n";
						TeukolskyMode teuk(-2, l, m, 0, n, geo);
						teuk.generateSolutions(geo);
						Zlmn = teuk.getTeukolskyAmplitude(Up);
						file << std::scientific << std::setprecision(15);
						file << p << "\t" << e << "\t" << std::real(Zlmn) <<"\t" << std::imag(Zlmn) << "\n";
					}
				}
				file.close();
			}
		}
	}

}

void test_teukolsky_mode(int s, int L, int m, int k, int n, GeodesicSource geo){
	TeukolskyMode teuk(s, L, m, k, n, geo);
	teuk.generateSolutions(geo);
	ComplexVector Rup = teuk.getHomogeneousRadialSolution(Up);
	ComplexVector Rin = teuk.getHomogeneousRadialSolution(In);
	Vector Sjm = teuk.getPolarSolution();
	std::cout << "ZlmUp = " << teuk.getTeukolskyAmplitude(Up) << "\n";
	std::cout << "ZlmIn = " << teuk.getTeukolskyAmplitude(In) << "\n";
	std::cout << "freq = " << teuk.getFrequency() << " \n";
	std::cout << "Rup(r = "<< teuk.getRadialPoints()[0] <<") = " << Rup[0]  << "\n";
	std::cout << "Rin(r = "<< teuk.getRadialPoints()[Rup.size()-1] <<") = " << Rin[Rup.size()-1]  << "\n";
	std::cout << "Sjm(theta = "<< teuk.getPolarPoints()[Sjm.size()-1] <<") = " << Sjm[Sjm.size()-1]  << "\n";
}

void test_teukolsky_data_generation(const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		std::cout << "Creating directory " << dir << " in "<< boost::filesystem::current_path() << "\n";
		boost::filesystem::create_directory(dir);
	}

	int rSampleNum = 10001;
	int omegaSampleNum = 5001;

	double a = 0.9;
	double separatrix = 2.320883041761887;
	double rmin = separatrix + 0.01;
	double rmax = 100.;
	double deltaU = 1./(rSampleNum - 1.); // for log scaling

	Vector rPts(rSampleNum);
	for(int i = 0; i < rSampleNum; i++){
		// rPts[i] = rmin + i*deltaR;
		rPts[i] = (pow(10., deltaU*i)*(rmax - rmin) + 10.*rmin - rmax)/9.; // log spacing
	}

	double omegaMin = 0.005;
	double omegaMax = 10. + omegaMin;
	double deltaOmega = (omegaMax - omegaMin)/(omegaSampleNum - 1.);

	Vector omegaPts(omegaSampleNum);
	for(int i = 0; i < omegaSampleNum; i++){
		omegaPts[i] = omegaMin + i*deltaOmega;
	}

	int s = 0;
	int l = 10;
	int m = 2;

	std::cout << "===============================\n";
	std::cout << "OMEGA = " << 1.e-4 << "\n";
	std::cout << "===============================\n";
	test_teukolsky_data_mode_generation_2(a, s, l, m, 1.e-4, rPts, dir);
	for(int i = 0; i < omegaSampleNum; i++){
		std::cout << "===============================\n";
		std::cout << "OMEGA = " << omegaPts[i] << "\n";
		std::cout << "===============================\n";
		test_teukolsky_data_mode_generation_2(a, s, l, m, omegaPts[i], rPts, dir);
	}
}

void test_teukolsky_data_mode_generation(double a, int s, int l, int m, double omega, Vector rPts, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	char buff[7];
	std::snprintf(buff, 7, "%f", omega);
	std::string filename = "teuk_data_";
	filename += buff;
	filename += ".txt";
	std::string filepath = dir + "/" + filename;

	if(!boost::filesystem::exists(filepath)){
		std::ofstream file;
		file.open(filepath);

		file << "a\ts\tl\tm\tomega\n";
		file << a << "\t" << s << "\t" << l << "\t" << m << "\t" << omega << "\n\n";
		file << "r\tRe[Rin]\tIm[Rin]\tRe[Rin']\tIm[Rin']\tRe[Rup]\tIm[Rup]\tRe[Rup']\tIm[Rup']\n";
		RadialTeukolsky teukSolution(a, s, l, m, omega, rPts);
		// teukSolution.generateRetardedBoundaryConditions();
		teukSolution.generateSolutions(In, AUTO);
		teukSolution.generateSolutions(Up, GSN);

		ComplexVector Rin = teukSolution.getSolution(In);
		ComplexVector RinP = teukSolution.getDerivative(In);
		ComplexVector Rup = teukSolution.getSolution(Up);
		ComplexVector RupP = teukSolution.getDerivative(Up);

		for(size_t i = 0; i < rPts.size(); i++){
			file << std::scientific << std::setprecision(15);
			file << rPts[i] << "\t" << std::real(Rin[i]) << "\t" << std::imag(Rin[i]) <<
				"\t" << std::real(RinP[i]) << "\t" << std::imag(RinP[i]) <<
					"\t" << std::real(Rup[i]) << "\t" << std::imag(Rup[i]) <<
						"\t" << std::real(RupP[i]) << "\t" << std::imag(RupP[i]) << "\n";
		}

		file.close();
	}
}

void test_teukolsky_data_mode_generation_2(double a, int s, int l, int m, double omega, Vector rPts, const std::string& dir){
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	char buff[7];
	std::snprintf(buff, 7, "%f", omega);
	std::string filename = "teuk_data_";
	filename += buff;
	filename += ".txt";
	std::string filepath = dir + "/" + filename;

	if(!boost::filesystem::exists(filepath)){
		std::ofstream file;
		file.open(filepath);

		file << "a\ts\tl\tm\tomega\n";
		file << a << "\t" << s << "\t" << l << "\t" << m << "\t" << omega << "\n\n";
		file << "r\tAbs[Rin]\tArg[Rin]\tAbs[Rin']\tArg[Rin']\tAbs[Rup]\tArg[Rup]\tAbs[Rup']\tArg[Rup']\n";
		RadialTeukolsky teukSolution(a, s, l, m, omega, rPts);
		// teukSolution.generateRetardedBoundaryConditions();
		if(s == 0){
			teukSolution.generateSolutions(AUTO);
		}else if(s < 0){
			teukSolution.generateSolutions(In, AUTO);
			teukSolution.generateSolutions(Up, GSN);
		}else{
			teukSolution.generateSolutions(In, GSN);
			teukSolution.generateSolutions(Up, AUTO);
		}

		ComplexVector Rin = teukSolution.getSolution(In);
		ComplexVector RinP = teukSolution.getDerivative(In);
		ComplexVector Rup = teukSolution.getSolution(Up);
		ComplexVector RupP = teukSolution.getDerivative(Up);

		for(size_t i = 0; i < rPts.size(); i++){
			file << std::scientific << std::setprecision(15);
			file << rPts[i] << "\t" << abs(Rin[i]) << "\t" << arg(Rin[i]) <<
				"\t" << abs(RinP[i]) << "\t" << arg(RinP[i]) <<
					"\t" << abs(Rup[i]) << "\t" << arg(Rup[i]) <<
						"\t" << abs(RupP[i]) << "\t" << arg(RupP[i]) << "\n";
		}

		file.close();
	}
}

void test_adiabatic_data_generation(){
	Vector rPts(201);
	for(int i = 0; i < 201; i++){
		rPts[i] = 3.5 + i/20.;
	}

	generate_adiabatic_circ_data(0.9, -1, 20, "inspiral");

	for(int i = 1; i < 9; i++){
		generate_adiabatic_circ_data(0.1*i, -1, 20, "inspiral");
	}

	generate_adiabatic_circ_data(0, -1, 20, "inspiral");

	for(int i = 1; i < 5; i++){
		generate_adiabatic_circ_data(0.9 + 0.02*i, -1, 20, "inspiral");
	}

	generate_adiabatic_circ_data(0.99, -1, 20, "inspiral");
	generate_adiabatic_circ_data(0.995, -1, 20, "inspiral");
	generate_adiabatic_circ_data(0.999, -1, 20, "inspiral");

}

void generate_adiabatic_inspiral_data(double a, int parallelFlag){
	int lmax = 15;
	int sampleSize_radius = 201;
	std::string dir = "../GEORG-Data/inspiral_temp/";
	if(!boost::filesystem::exists(dir)){
		boost::filesystem::create_directory(dir);
	}

	if(a >= 0){
		dir += "a" + std::to_string(int(abs(a)*10000));
	}else{
		dir += "an" + std::to_string(int(abs(a)*10000));
	}
	if(parallelFlag == 1){
		output_adiabatic_circ_spherical_mode_data_radial_sample_parallel(a, lmax, sampleSize_radius, dir);
	}else{
		output_adiabatic_circ_spherical_mode_data_radial_sample(a, lmax, sampleSize_radius, dir);
	}
	generate_adiabatic_circ_remainder_data_radial_sample(a, lmax + 1, sampleSize_radius, dir);
	generate_adiabatic_circ_trajectory(dir);
}
