from kerrgeocy import KerrGeodesicWrapper
import kerrgeocy
import numpy as np

kerr_orbital_constants = kerrgeocy.kerr_orbital_constants_wrapper
kerr_mino_frequencies = kerrgeocy.kerr_mino_frequencies_wrapper

def kerr_fundamental_frequencies(a, p, e, x):
    ups = kerr_mino_frequencies(a, p, e, x)
    return ups[1:]/ups[0]

class KerrGeodesic(KerrGeodesicWrapper):
    def __init__(self, a, p, e, x, nsamples=2**8):
        KerrGeodesicWrapper.__init__(self, a, p, e, x, nsamples=nsamples)

        # unfortunately we are going to waste some memory copying over these vectors so that
        # they are quickly accessible as numpy arrays. Calling time_accumulation and casting
        # C++ vectors to numpy arrays takes {\mu}s while just calling numpy arrays takes ns
        time_accumulation_r = self.get_time_accumulation(1)
        time_accumulation_th = self.get_time_accumulation(2)
        self._time_accumulation = {
            1: time_accumulation_r,
            2: time_accumulation_th
        }
        phi_accumulation_r = self.get_azimuthal_accumulation(1)
        phi_accumulation_th = self.get_azimuthal_accumulation(2)
        self._azimuthal_accumulation = {
            1: phi_accumulation_r,
            2: phi_accumulation_th
        }

        rp = self.get_radial_points()
        thp = self.get_polar_points()
        self._radial_points = rp
        self._polar_points = thp

    @property
    def time_accumulation(self):
        return self._time_accumulation

    @property
    def radial_points(self):
        return self._radial_points

    @property
    def polar_points(self):
        return self._polar_points
    
    @property
    def azimuthal_accumulation(self):
        return self._azimuthal_accumulation

    def __call__(self, la):
        if isinstance(la, np.ndarray):
            return self.position_vec(la)
        else:
            return self.position(la)
            