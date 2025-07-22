








# ****************This file has all the function used to compute fluxes from Pybhpt code in parallel 

import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from pybhpt.hertz import HertzMode
from pybhpt.hertz import available_gauges
from pybhpt.flux import FluxMode, FluxList
from concurrent.futures import ProcessPoolExecutor
import os
from few.utils.geodesic import get_separatrix
from numpy import sin, cos,pi



def pow(x, y=2):
    return np.power(x, y)
def sqrt(x):
    return np.sqrt(x)
def exp(x):
    return np.exp(x)
def ISCO(a):
    term1 = 3 + sqrt(3 * pow(a) + pow(1 + pow(1 - pow(a), 1/3) * (pow(1 - a, 1/3) + pow(1 + a, 1/3))))
    term2 = sqrt((2 - pow(1 - pow(a), 1/3) * (pow(1 - a, 1/3) + pow(1 + a, 1/3))) * 
                      (4 + pow(1 - pow(a), 1/3) * (pow(1 - a, 1/3) + pow(1 + a, 1/3)) + 
                       2 * sqrt(3 * pow(a) + pow(1 + pow(1 - pow(a), 1/3) * (pow(1 - a, 1/3) + pow(1 + a, 1/3))))))

    return term1 - term2





def Kerr_circ_flux_calc(a, p, lmax):
    # if p<3.5:
    #     lmax = 30
    # elif (p>3.5)&(p<8.0):
    #     lmax = 20
    # else: lmax = 15
#     lmax = 20        
    e = 0.0
    x = 1.0
    if a<0.0:
        a = -a
        x = -1.0
    nsamples = 2**8
    geo = KerrGeodesic(a, p, e, x, nsamples)
    omega_phi = geo.frequencies[2]
    Ehor  = np.zeros(1)
    Einf = np.zeros(1)
    Lhor  = np.zeros(1)
    Linf = np.zeros(1)
    s = -2
    k = 0
    n = 0
    print(lmax)
    for l in range(abs(s),lmax+1):
        for m in range(1,l+1):
            start_time = time.time()
            teuk = TeukolskyMode(-2, l, m, k, n, geo)
            teuk.solve(geo)
            Ehor += 2*FluxMode(geo, teuk).horizonfluxes[0]
            Einf += 2*FluxMode(geo, teuk).infinityfluxes[0]
            Lhor += 2*FluxMode(geo, teuk).horizonfluxes[1]
            Linf += 2*FluxMode(geo, teuk).infinityfluxes[1]
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"Elapsed Time: {elapsed_time} seconds")
    Edot_tot = Einf[0] + Ehor[0]
    Ldot_tot = Linf[0] + Lhor[0]         
    return a, p, omega_phi, Einf[0], Ehor[0], Edot_tot, Linf[0], Lhor[0], Ldot_tot


















def compute_terms(args):
    a, p, l, m, k, n, = args
    nsamples = 2**4 # Again I am running it with lower nsamples #In the new run I am hanging it to 2^9 , last run was 2^8
    ########### For circular orbits, I dont need many nsampels, so I am using 2^4.tested alady compared to 2**9 no difference.
    e = 0.0
    x = 1.0
    if a<0.0:
        a = -a
        x = -1.0
    geo = KerrGeodesic(a, p, e, x, nsamples)
    omega_phi = geo.frequencies[2]
    teuk = TeukolskyMode(-2, l, m, k, n, geo)
    # teuk.solve(geo,method = "AUTO", nsamples = nsamples)
    Ehor = 2*FluxMode(geo, teuk).horizonfluxes[0]
    Einf = 2*FluxMode(geo, teuk).infinityfluxes[0]
    Lhor = 2*FluxMode(geo, teuk).horizonfluxes[1]
    Linf = 2*FluxMode(geo, teuk).infinityfluxes[1]
    return Ehor, Einf, Lhor, Linf, omega_phi, l, m

def Kerr_circ_flux_calc_parallel(a, p, lmax, save_modes_option = False, path=None,  grid_option_u = "normal"):
    # x = 1.0
    # if a<0.0:
    #     x = -1.0
    x = 1.0 #np.sign(a)    # in the new version of FEW the sign of a is not needed for x, a<0 is automatically taken care of
    ps = get_separatrix(a,0.0,x)#ISCO(a) ps takes the abs value of spin, so I don't need to do a = -a in case of retrograde
    u = np.log(p - ps + 3.9)
    Ehor  = np.zeros(1)
    Einf = np.zeros(1)
    Lhor  = np.zeros(1)
    Linf = np.zeros(1)
    s = -2
    k = 0
    n = 0
    tasks = [(a,p,l, m, k, n) for l in range(abs(s), lmax+1) for m in range(1, l+1)]
    num_tasks = len(tasks)
    modes_array = np.zeros((num_tasks,12))
    # num_cores = 40
    avail_cores = os.cpu_count()
    num_cores = int(avail_cores) #40
    # print("number of cores:", num_cores, num_tasks)
    with ProcessPoolExecutor(max_workers = num_cores) as executor:
        for idx,result in enumerate(executor.map(compute_terms, tasks)):
            Ehor_lm, Einf_lm, Lhor_lm, Linf_lm, omega_phi, l, m = result
            Edot_tot_lm = Einf_lm + Ehor_lm
            Ldot_tot_lm = Linf_lm + Lhor_lm
            Ehor += result[0]
            Einf += result[1]
            Lhor += result[2]
            Linf += result[3]
            omega_phi = result[4]
            # print(idx)
            modes_array[idx] = np.array([l, m,a, p, u, omega_phi, Einf_lm, Ehor_lm, Edot_tot_lm, Linf_lm, Lhor_lm, Ldot_tot_lm])
    Edot_tot = Einf[0] + Ehor[0]
    Ldot_tot = Linf[0] + Lhor[0]
    if(grid_option_u == "Cheby"):
        u = p - ps + 1.
    # u = p - ps + 3.9 #This one is for the Cheby grid that I am choosing
    if(save_modes_option):
        if path is None:
            print("Please provide a path to save the modes")

        dir_path = path  + "/modes/"
        os.makedirs(dir_path, exist_ok=True)
        a_formatted = f"{a:.5f}"
        file_path = os.path.join(dir_path, f"modes_array_{a:.6f}_{p}.npy")
        np.save(file_path, modes_array)
        # stores the modes array in a file in the same directory

    return a, p, u, omega_phi, Einf[0], Ehor[0], Edot_tot, Linf[0], Lhor[0], Ldot_tot


