// utils.hpp

/* A very basis header file that includes the standard C++ headers so that they
do not need to be individually added to all the other user-defined header files
*/

#include "utils.hpp"

Vector downsample_base_2(Vector v, int sampleRate){
  int size;
  if(v.size() % 2 == 1){
    size = (v.size() - 1)/sampleRate + 1;
  }else{
    size = v.size()/sampleRate;
  }
  Vector w(size);
  for(size_t i = 0; i < w.size(); i++){
    w[i] = v[i*sampleRate];
  }

  return w;
}

Vector unwrap_phase(const Vector& phase){
	Vector phase_unwrap(phase.size());
	phase_unwrap[0] = phase[0];
	for(size_t i = 1; i < phase_unwrap.size(); i++){
		double delta_phase = phase[i] - phase[i-1];
		if(delta_phase > M_PI){
			delta_phase -= 2.*M_PI;
		}else if(delta_phase < - M_PI){
			delta_phase += 2.*M_PI;
		}
		phase_unwrap[i] = phase_unwrap[i - 1] + delta_phase;
	}

	return phase_unwrap;
}

void cubic_solver(double &rp, double &ra, double &r3, double A3, double A2, double A1, double A0) {
    double a_c = 1.0; 
    double b = A2 / A3;
    double c = A1 / A3;
    double d = A0 / A3;

    double p_c = (3.0 * c - b * b) / 3.0;
    double q_c = (2.0 * std::pow(b, 3) - 9.0 * b * c + 27.0 * d) / 27.0;
    double disc = (q_c * q_c / 4.0) + (std::pow(p_c, 3) / 27.0);

    if (disc >= 0) {
        double sqrt_disc = std::sqrt(disc);
        double u = std::cbrt(-q_c / 2.0 + sqrt_disc);
        double v = std::cbrt(-q_c / 2.0 - sqrt_disc);
        r3 = u + v - b / 3.0;
        ra = -0.5 * (u + v) - b / 3.0;
        rp = ra; 
    } else {
        double r = std::sqrt(-p_c / 3.0);
        double theta = std::acos(-q_c / (2.0 * r * r * r));
        ra = 2.0 * r * std::cos(theta / 3.0) - b / 3.0;
        r3 = 2.0 * r * std::cos((theta + 2.0 * M_PI) / 3.0) - b / 3.0;
        rp = 2.0 * r * std::cos((theta - 2.0 * M_PI) / 3.0) - b / 3.0;
    }
}
