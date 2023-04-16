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
