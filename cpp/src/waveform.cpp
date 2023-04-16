// waveform.cpp

#include "waveform.hpp"

#define WAVEFORM_EPSILON 1.e-6
#define KN_MAX 200

WaveformModeLMKN generateWaveformModeLMKN(int, int m, int k, int n, GeodesicSource source, TeukolskyMode teuk){
	Complex amplitude;
	double frequency = source.getTimeFrequency(m, k, n);
	if(source.getEccentricity() == 0. && abs(n) > 0){
		amplitude = 0.;
	}else if(source.getInclination() == 0. && abs(k) > 0){
		amplitude = 0.;
	}else if(frequency == 0.){
		amplitude = 0.;
	}else{
		teuk.generateSolutions(source);
		amplitude = -teuk.getTeukolskyAmplitude(Up)*pow(teuk.getFrequency(), -2);
	}
	WaveformModeLMKN waveform = {
		.k = k,
		.n = n,
		.amplitude = amplitude,
		.frequency = frequency};
	return waveform;
}

WaveformModeLMKN generateWaveformModeLMKN(int L, int m, int k, int n, GeodesicSource source){
	Complex amplitude;
	double frequency = source.getTimeFrequency(m, k, n);
	if(source.getEccentricity() == 0. && abs(n) > 0){
		amplitude = 0.;
	}else if(source.getInclination() == 0. && abs(k) > 0){
		amplitude = 0.;
	}else if(frequency == 0.){
		amplitude = 0.;
	}else{
		TeukolskyMode teuk(L, m, k, n, source);
		teuk.generateSolutions(source);
		amplitude = -teuk.getTeukolskyAmplitude(Up)*pow(teuk.getFrequency(), -2);
	}
	WaveformModeLMKN waveform = {
		.k = k,
		.n = n,
		.amplitude = amplitude,
		.frequency = frequency};
	return waveform;
}

WaveformModeLMKN generateWaveformModeLMKN(int L, int m, int k, int n, GeodesicSource source, TeukolskyMode teuk, double theta){
	Complex amplitude;
	double frequency = source.getTimeFrequency(m, k, n);
	if(source.getEccentricity() == 0. && abs(n) > 0){
		amplitude = 0.;
	}else if(source.getInclination() == 0. && abs(k) > 0){
		amplitude = 0.;
	}else if(frequency == 0.){
		amplitude = 0.;
	}else{
		Vector thp(1, theta);
		SpinWeightedHarmonic swsh(-2, L, m, source.getBlackHoleSpin()*frequency, thp);

		teuk.generateSolutions(source);
		amplitude = -teuk.getTeukolskyAmplitude(Up)*pow(teuk.getFrequency(), -2)*swsh.getSolution()[0];
	}
	WaveformModeLMKN waveform = {
		.k = k,
		.n = n,
		.amplitude = amplitude,
		.frequency = frequency};
	return waveform;
}

WaveformModeLMKN generateWaveformModeLMKN(int L, int m, int k, int n, GeodesicSource source, double theta){
	Complex amplitude;
	double frequency = source.getTimeFrequency(m, k, n);
	if(source.getEccentricity() == 0. && abs(n) > 0){
		amplitude = 0.;
	}else if(source.getInclination() == 0. && abs(k) > 0){
		amplitude = 0.;
	}else if(frequency == 0.){
		amplitude = 0.;
	}else{
		Vector thp(1, theta);
		SpinWeightedHarmonic swsh(-2, L, m, source.getBlackHoleSpin()*frequency, thp);

		TeukolskyMode teuk(L, m, k, n, source);
		teuk.generateSolutions(source);
		amplitude = -teuk.getTeukolskyAmplitude(Up)*pow(teuk.getFrequency(), -2)*swsh.getSolution()[0];
	}
	WaveformModeLMKN waveform = {
		.k = k,
		.n = n,
		.amplitude = amplitude,
		.frequency = frequency};
	return waveform;
}

WaveformModeLM::WaveformModeLM(): _L(2), _m(2), _nmax(0), _kmax(0), _kmin(0), _lmknList(KN_MAX){}
WaveformModeLM::WaveformModeLM(int L, int m): _L(L), _m(m), _nmax(0), _kmax(0), _kmin(0), _lmknList(KN_MAX){}
WaveformModeLM::~WaveformModeLM(){}

int WaveformModeLM::generateModes(GeodesicSource source, double theta){
	int globalKMin = _kmin, globalKMax = _kmax;
	int kn = 0;
	WaveformModeLMKN waveformMode = generateWaveformModeLMKN(_L, _m, _kmax, _nmax, source, theta);
	_lmknList[kn] = waveformMode;
	_power = 0.;
	double modePower = abs(waveformMode.amplitude);
	_power += modePower;
	if(source.getEccentricity() == 0. && abs(source.getInclination()) == 1.){
		_kmax = globalKMax;
		_kmin = globalKMin;
		_lmknList.resize(kn + 1);
		//std::cout << "WAVEFORM: Generated all necessary modes between " << _kmin << " <= k <= " << _kmax << " and " << _nmin << " <= n <= " << _nmax << "\n";
		return 0;
	}

	while(modePower/_power > WAVEFORM_EPSILON && kn < KN_MAX){
		_kmax = 0;
		_kmin = 0;
		while(modePower/_power > WAVEFORM_EPSILON && kn < KN_MAX){
			_kmax++;
			globalKMax = globalKMax > _kmax? globalKMax : _kmax;
			kn++;
			waveformMode = generateWaveformModeLMKN(_L, _m, _kmax, _nmax, source, theta);
			_lmknList[kn] = waveformMode;
			modePower = abs(waveformMode.amplitude);
			_power += modePower;
		}

		_kmin --;
		globalKMin = globalKMin < _kmin? globalKMin : _kmin;
		kn++;
		waveformMode = generateWaveformModeLMKN(_L, _m, _kmax, _nmax, source, theta);
		_lmknList[kn] = waveformMode;
		modePower = abs(waveformMode.amplitude);
		_power += modePower;
		while(modePower/_power > WAVEFORM_EPSILON && kn < KN_MAX){
			_kmin--;
			globalKMin = globalKMin < _kmin? globalKMin : _kmin;
			kn++;
			waveformMode = generateWaveformModeLMKN(_L, _m, _kmax, _nmax, source, theta);
			_lmknList[kn] = waveformMode;
			modePower = abs(waveformMode.amplitude);
			_power += modePower;
		}
		_nmax++;
	}

	_kmax = globalKMax;
	_kmin = globalKMin;
	_lmknList.resize(kn + 1);

	return 0;
}

double WaveformModeLM::getPower(){ return _power; }

void WaveformModeLM::printPowerSpectrum(){
	std::cout << "n/k \t";
	int kn = 0;
	Vector kModeZero(_kmax - _kmin + 1, 0.);
	Vector kMode = kModeZero;
	for(int k = _kmin; k < _kmax; k++){
		std::cout << k << "\t";
	}
	std::cout << _kmax << "\n";
	for(int n = 0; n <= _nmax; n++){
		std::cout << n << "\t";
		for(int k = 0; k <= _kmax; k++){
			if(_lmknList[kn].k == k){
				kMode[k] = abs(_lmknList[kn].amplitude);
				kn++;
			}
		}

		for(int k = -1; k >= _kmin; k--){
			if(_lmknList[kn].k == k){
				kMode[_kmax - k] = abs(_lmknList[kn].amplitude);
				kn++;
			}
		}

		for(int k = _kmin; k < 0; k++){
			std::cout << kMode[_kmax - k] << "\t";
		}

		for(int k = 0; k < _kmax; k++){
			std::cout << kMode[k] << "\t";
		}

		std::cout << kMode[_kmax] << "\n";
		kMode = kModeZero;
	}
}

ComplexVector WaveformModeLM::getWaveform(Vector t){
	ComplexVector waveform(t.size());
	for(size_t i = 0; i < t.size(); i++){
		waveform[i] = 0.;
		for(size_t j = 0; j < _lmknList.size(); j++){
			waveform[i] += _lmknList[j].amplitude*exp(-I*_lmknList[j].frequency*t[i]);
		}
	}

	return waveform;
}

void WaveformModeLM::writeWaveform(Vector t, const std::string& fileName){
	std::cout << "WAVEFORM: Writing waveform to file " << fileName << "\n";
	std::ofstream file;
	file.open(fileName);
	file << "t" << "\t" << "Re[h]" << "\t" << "Im[h]" << "\n";
	file << std::setprecision(6);
	ComplexVector waveform = getWaveform(t);
	for(size_t i = 0; i < waveform.size(); i++){
		file << t[i] << "\t" << std::real(waveform[i]) << "\t" << std::imag(waveform[i]) << "\n";
	}

	file.close();
}

Waveform::Waveform(int Lmax): _Lmax(Lmax), _harmonicPower(_Lmax - 2, 0.), _lmList((Lmax + 3)*(Lmax - 1)) {}
Waveform::~Waveform(){}

int Waveform::generateModes(GeodesicSource source, double theta){
	_a = source.getBlackHoleSpin();
	_observerAngle = theta;
	int lmCounter = 0;
	int l = 2;
	double power = 0., modePower = 0.;
	WaveformModeLM waveformMode;
	for(int m = -_Lmax; m <= _Lmax; m++){
		if(abs(m) < 2) l = 2; else l = abs(m);
		waveformMode = WaveformModeLM(l, m);
		waveformMode.generateModes(source, theta);
		_lmList[lmCounter] = waveformMode;
	    modePower = waveformMode.getPower();
	   	power += modePower;
		_harmonicPower[l - 2] += modePower;
		lmCounter++;
		l++;
		while(modePower/power > WAVEFORM_EPSILON && l <= _Lmax){
			waveformMode = WaveformModeLM(l, m);
			waveformMode.generateModes(source, theta);
			_lmList[lmCounter] = waveformMode;
		    modePower = waveformMode.getPower();
		   	power += modePower;
			_harmonicPower[l - 2] += modePower;
			lmCounter++;
			l++;
		}
	}
	_lmList.resize(lmCounter);
	std::cout << "Relative power of l = " << _Lmax << " mode is " << _harmonicPower[_Lmax - 2]/power << "\n";

	return 0;
}

ComplexVector Waveform::getWaveform(Vector t){
	ComplexVector waveform(t.size(), 0.);
	ComplexVector waveformMode(t.size());
	for(size_t lm = 0; lm < _lmList.size(); lm++){
		waveformMode = _lmList[lm].getWaveform(t);
		for(size_t i = 0; i < t.size(); i++){
			waveform[i] += waveformMode[i];
		}
	}

	return waveform;
}

Vector Waveform::getHarmonicPower(){ return _harmonicPower; }
double Waveform::getHarmonicPower(int L){
	if(L > _Lmax){
		std::cout << "Error: Modes only generated up to l = " << _Lmax << "\n";
		return 0.;
	}

	return _harmonicPower[L - 2];
}

void Waveform::writeWaveform(Vector t, const std::string& fileName){
	std::cout << "WAVEFORM: Writing waveform to file " << fileName << "\n";
	std::ofstream file;
	file.open(fileName);
	file << "t" << "\t" << "Re[h]" << "\t" << "Im[h]" << "\n";
	file << std::setprecision(6);
	ComplexVector waveform = getWaveform(t);
	for(size_t i = 0; i < waveform.size(); i++){
		file << t[i] << "\t" << std::real(waveform[i]) << "\t" << std::imag(waveform[i]) << "\n";
	}

	file.close();
}
