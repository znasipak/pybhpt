
import numpy as np
import matplotlib.pyplot as plt

from pybhpt.swsh import SWSH
from pybhpt.geo import KerrGeodesic
from pybhpt.teuk import TeukolskyMode
from pybhpt.hertz import HertzMode
from pybhpt.hertz import available_gauges
from pybhpt.flux import FluxMode, FluxList
import signal
import sys


from concurrent.futures import ProcessPoolExecutor
import os


# def signal_handler(sig, frame):
#     print('Interrupt received, terminating pool...')
#     executor.shutdown(wait=False, cancel_futures=True)
#     sys.exit(0)

# def increase_geo_sampling(e, p, nsamples_updated):
#     a = 0.0
#     x = 1.0
#     return KerrGeodesic(a, p, e, x, nsamples=nsamples_updated)

def init_worker(geo_args):
    global GEO_LOW, GEO_MED, GEO_HIGH, GEO_HIGH2, GEO_HIGH3, GEO_HIGH4, GEO_HIGH5, GEO_HIGH6
    a, p, e, x, base_nsamples = geo_args
    GEO_LOW = KerrGeodesic(a, p, e, x, base_nsamples)
    GEO_MED = KerrGeodesic(a, p, e, x, base_nsamples*2)  # 2x resolution
    GEO_HIGH = KerrGeodesic(a, p, e, x, base_nsamples*2**2)  # 4x resolution
    GEO_HIGH2 = KerrGeodesic(a, p, e, x, base_nsamples*2**3)  # 8x resolution
    GEO_HIGH3 = KerrGeodesic(a, p, e, x, base_nsamples*2**4)  # 16x resolution
    GEO_HIGH4 = KerrGeodesic(a, p, e, x, base_nsamples*2**5)  # 32x resolution
    GEO_HIGH5 = KerrGeodesic(a, p, e, x, base_nsamples*2**6)  # 64x resolution
    GEO_HIGH6 = KerrGeodesic(a, p, e, x, base_nsamples*2**7)  # 128x resolution
    # global GEO
    # GEO = geo_arg



def compute_terms(args):
    e, p, l, m, k, n,_ = args
    # nsamples = 2**4 # Again I am running it with lower nsamples #In the new run I am hanging it to 2^9 , last run was 2^8
    a = 0.0
    x = 1.0


    if n < 0:
        geo = GEO_HIGH5
    elif n < 10:
        geo = GEO_LOW
    elif n < 30:
        geo = GEO_MED
    elif n < 60:
        geo = GEO_HIGH
    elif n < 100:
        geo = GEO_HIGH2
    elif abs(n) < 150:
        geo = GEO_HIGH3
    elif abs(n) < 200:
        geo = GEO_HIGH4
    elif abs(n) < 300:
        geo = GEO_HIGH5
    elif abs(n) < 400:
        geo = GEO_HIGH6
    else:
        geo = GEO_HIGH6            


    
    # if abs(n) > 84 and ns < 2**11:
    #     ns = 2**11
    #     geo = KerrGeodesic(a, p, e, x, ns)
  
    # print("nsamples:", ns)


    omega_phi = geo.frequencies[2]
    teuk = TeukolskyMode(-2, l, m, k, n, geo)
    # teuk.solve(geo,method = "AUTO", nsamples = nsamples)
    Ehor = 2*FluxMode(geo, teuk).horizonfluxes[0]
    Einf = 2*FluxMode(geo, teuk).infinityfluxes[0]
    Lhor = 2*FluxMode(geo, teuk).horizonfluxes[1]
    Linf = 2*FluxMode(geo, teuk).infinityfluxes[1]
    return Ehor, Einf, Lhor, Linf, omega_phi, l, m, n

def Schw_Ecc_flux_calc_parallel(e, p, lmax, nr_min, nr_max, nsamples, save_modes_option = False, path=None,  grid_option_u = "normal"):
    # x = 1.0
    # if a<0.0:
    #     x = -1.0
    x = 1.0 #np.sign(a)    # in the new version of FEW the sign of a is not needed for x, a<0 is automatically taken care of
    # ps = get_separatrix(a,0.0,x)#ISCO(a) ps takes the abs value of spin, so I don't need to do a = -a in case of retrograde
    # u = np.log(p - ps + 3.9)
    Ehor  = np.zeros(1)
    Einf = np.zeros(1)
    Lhor  = np.zeros(1)
    Linf = np.zeros(1)
    s = -2
    k = 0
    a = 0.0
    rp = p / (1.0 + e) 
    # geo = KerrGeodesic(a, p, e, x, nsamples)
    geo_args = (a, p, e, x, nsamples)



    if(save_modes_option):
        if path is None:
            print("Please provide a path to save the modes")
        dir_path = path  + "/Results/modes_grid_full/"
        os.makedirs(dir_path, exist_ok=True)
        file_path = os.path.join(dir_path, f"modes_array_{e:.4f}_{rp:.1f}.txt")
        fmt = "%d %d %d %.5f %.5f %.5f %1.12e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e"
        with open(file_path, 'a') as f:
            f.write("# l m n e rp p omega_phi Einf Ehor Edot_tot Linf Lhor Ldot_tot\n")
        
    
    # n = 0
    # tasks = [(e,p,l, m, k, n, nsamples) for l in range(abs(s), lmax+1) for m in range(0, l+1) for n in range(nr_min, nr_max+1) ]
    # tasks = [(e,p,l, m, k, n, nsamples)  for m in range(5, lmax+1) for l in range(max(abs(s), m), lmax+1) for n in range(nr_min, nr_max+1) ]
    tasks = [(e,p,l, m, k, n, nsamples)  for (l,m) in zip(range(2, lmax+1), range(2, lmax+1)) for n in range(nr_min, nr_max+1) ]  #### onl l=m modes
    # tasks = [(e,p,l, m, k, n, nsamples)  for (l,m) in zip(range(10, lmax+1), range(10, lmax+1)) for n in range(nr_min, nr_max+1) ]
    num_tasks = len(tasks)
    # modes_array = np.zeros((num_tasks,13))
    # num_cores = 40
    avail_cores = 60#os.cpu_count()
    num_cores = int(avail_cores) #40
    print("number of cores:", num_cores, num_tasks)

    chunk_size = 100
    write_buffer  = np.zeros((chunk_size, 13))
    buffer_idx = 0
    row_template = np.zeros(13)
    with open(file_path, 'a') as f:
        with ProcessPoolExecutor(max_workers = num_cores, initializer=init_worker, initargs=(geo_args,)) as executor:
            for idx,result in enumerate(executor.map(compute_terms, tasks)):
                Ehor_lm, Einf_lm, Lhor_lm, Linf_lm, omega_phi, l, m, n = result
                Edot_tot_lm = Einf_lm + Ehor_lm
                Ldot_tot_lm = Linf_lm + Lhor_lm
                Ehor += result[0]
                Einf += result[1]
                Lhor += result[2]
                Linf += result[3]
                omega_phi = result[4]
                # modes_array[idx] = np.array([l, m, n, e, rp, p, omega_phi, Einf_lm, Ehor_lm, Edot_tot_lm, Linf_lm, Lhor_lm, Ldot_tot_lm])
                # np.savetxt(f, modes_array[idx][None, :], fmt=fmt)
                row_template[:] = [l, m, n, e, rp, p, omega_phi, Einf_lm, Ehor_lm, Edot_tot_lm, Linf_lm, Lhor_lm, Ldot_tot_lm]
                write_buffer[buffer_idx] = row_template
                buffer_idx += 1
                if buffer_idx == chunk_size:
                    np.savetxt(f, write_buffer, fmt=fmt)
                    f.flush()
                    os.fsync(f.fileno())
                    buffer_idx = 0  # reset the buffer index

        if buffer_idx > 0:
            np.savetxt(f, write_buffer[:buffer_idx], fmt=fmt)
            f.flush()
            os.fsync(f.fileno())


    Edot_tot = Einf[0] + Ehor[0]
    Ldot_tot = Linf[0] + Lhor[0]

    


    return e, rp, p, Einf[0], Ehor[0], Edot_tot, Linf[0], Lhor[0], Ldot_tot







import time


t_init = time.time()


nr_min = -80
nr_max = 800

if len(sys.argv) != 3:
    print("Usage: python your_script.py <rp> <ecc>")
    sys.exit(1)
rp = float(sys.argv[1])
ecc = float(sys.argv[2])


# rp = 6.0
# ecc = 0.94
p = rp * (1 + ecc)
save_modes_option = True
lmax = 10
current_path = "/home/hkhalvati/Downloads/Teuk_python_zach/High_Eccetric_Pybhpt/"
# current_path = "/mnt/beegfs/hkhalvati/High_Eccetric_Pybhpt/"

nsamples = 2**9


Schw_Ecc_flux_calc_parallel(ecc, p, lmax, nr_min, nr_max,nsamples, save_modes_option = save_modes_option, path=current_path)



t_final = time.time()
print("Time taken:", t_final - t_init, "seconds")
