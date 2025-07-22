
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys 
# for i in range(0,99,3):
#     print("echo -e","\"",i,"\\n",i+3,"\\n","\"","| python run_flux.py")


# path = '/home/hkhalvati/Downloads/Teuk_python_zach/'
current_path = os.getcwd()
path = current_path + '/..'
print(path)
sys.path.append(path)

from src.flux_pybhpt import pow, ISCO, exp, sqrt, Kerr_circ_flux_calc_parallel
from numpy import sin, cos,pi
from few.utils.geodesic import get_separatrix



ni_in = int(input("Enter the value for ni_in: "))
ni_f = int(input("Enter the value for ni_f: "))
# ni_in = int(input())
# ni_f = int(input())


fmt = '%2.12g' , '%2.12g', '%2.12g' , '%2.16g', '%2.16g' , '%2.16g', '%2.16g' , '%2.16g', '%2.16g' , '%2.16g'
# ni_f = 2
# ni_in = 1
chunck_size = ni_f - ni_in
lmax = 60


dir_path = f"../result_lmax{lmax}_log_both_fullgrid_bias_factor4/"
# Create the directory if it doesn't exist
os.makedirs(dir_path, exist_ok=True)


du = 0.025
nj = int((3.8186394 - 1.3686394)/du) + 1 #this goes up to larger distances for u
chunck_size = ni_f - ni_in
print(nj, chunck_size)

flux_data = np.zeros((chunck_size*nj,10))
new_result_array = np.zeros((chunck_size*nj,11))

######## Changing the logarithmic grid to be all from -0.99 to +0.99
n_alpha = 100
a1 = 0.99
a2 = -0.99
alpha1 = np.log(a1 + 2.5) 
alpha2 = np.log(a2 + 2.5) 

alpha = np.linspace(alpha1, alpha2, n_alpha)
bias_factor = 4  # Adjust for more density near a = 0.99
scaled_alpha = alpha1 + (alpha - alpha1) ** bias_factor / (alpha2 - alpha1) ** (bias_factor - 1)


#modified the format to include alpha0
fmt = '%2.12g' , '%2.12g', '%2.12g' , '%2.16g', '%2.16g' , '%2.16g', '%2.16g' , '%2.16g', '%2.16g' , '%2.16g', '%2.12g'




for i, alpha0 in enumerate(scaled_alpha[ni_in:ni_f]):
    a = np.exp(alpha0) - 2.5
    for j in range(nj):
        init_time = time.time()
        u = 1.3686394259 + du*j #<-----similar grid as mine-------#u = 1.37 + 0.025*j
        u = float("{:.14f}".format(u))
        # ps = ISCO(a) # ISCO function is not good for retro orbits
        x = 1.0#np.sign(a)     ##### in the new version of FEW the sign of a is not needed for x, a<0 is automatically taken care of
        # x = 1.0
        # if a <0.0:
        #     x = - 1.0
        ps = get_separatrix(a, 0.0,x)
        # breakpoint()
        p = exp(u) + ps - 3.9 
        flux_data[i*nj+j] = Kerr_circ_flux_calc_parallel(a, p, lmax )
        # breakpoint()
        new_result_array[i*nj+j] = np.hstack((flux_data[i*nj+j], alpha0))
        print(i, j, nj, i*nj+j, a,alpha0, p, ps, u)
        if (j)%5==0:
            np.savetxt(dir_path+"flux_pybhpt{}_{}_log_both.txt".format(ni_in, ni_f),new_result_array, fmt=fmt)
    np.savetxt(dir_path + "flux_pybhpt{}_{}_log_both.txt".format(ni_in, ni_f),new_result_array, fmt=fmt)
np.savetxt(dir_path + "flux_pybhpt{}_{}_log_both.txt".format(ni_in, ni_f),new_result_array, fmt=fmt)
