// waveform.hpp

#ifndef WAVEFORM_HPP
#define WAVEFORM_HPP

#include "teukolsky.hpp"

typedef struct WaveformModeLMKNStruct{
	int k;
	int n;
	Complex amplitude;
	double frequency;
} WaveformModeLMKN;

class WaveformModeLM{
public:
	WaveformModeLM();
	WaveformModeLM(int L, int m);
	~WaveformModeLM();

	int generateModes(GeodesicSource source, double theta);
	double getPower();
	void printPowerSpectrum();
	ComplexVector getWaveform(Vector t);

	void writeWaveform(Vector t, const std::string& input);

private:
	int _L;
	int _m;

	int _nmax;
	int _kmax;
	int _kmin;
	double _power;
	std::vector<WaveformModeLMKN> _lmknList;
};

class Waveform{
public:
	Waveform(int Lmax);
	~Waveform();

	int generateModes(GeodesicSource source, double theta);
	ComplexVector getWaveform(Vector t);
	ComplexVector getWaveform(int L, int m, Vector t);
	Vector getHarmonicPower();
	double getHarmonicPower(int L);

	void writeWaveform(Vector t, const std::string& input);
	void writeWaveform(int L, int m, Vector t, const std::string& input);

private:
	double _a;
	double _observerAngle;
	int _Lmax;

	Vector _harmonicPower;
	std::vector<WaveformModeLM> _lmList;
};

WaveformModeLMKN generateWaveformModeLMKN(int L, int m, int k, int n, GeodesicSource source);
WaveformModeLMKN generateWaveformModeLMKN(int L, int m, int k, int n, GeodesicSource source, TeukolskyMode teuk);

WaveformModeLMKN generateWaveformModeLMKN(int L, int m, int k, int n, GeodesicSource source, double theta);
WaveformModeLMKN generateWaveformModeLMKN(int L, int m, int k, int n, GeodesicSource source, TeukolskyMode teuk, double theta);

#endif
