import os
from Cython.Build import cythonize
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
from subprocess import run
import sys

compiler_flags = ["-std=c++11", "-fopenmp", "-march=native"]
libraries = ["gsl", "gslcblas", "boost_filesystem"]

# set some flags and libraries depending on the system platform
if sys.platform.startswith('win32'):
    compiler_flags.append('/Od')
    libraries.append('gomp')
elif sys.platform.startswith('darwin'):
    compiler_flags.append('-O2')
    libraries.append('omp')
elif sys.platform.startswith('linux'):
    compiler_flags.append('-O2')
    libraries.append('gomp')
    

CFLAGS = os.getenv("CFLAGS")
if CFLAGS is None:
    CFLAGS = ""
else:
    CFLAGS = str(CFLAGS)
os.environ["CFLAGS"] = CFLAGS
os.environ["CFLAGS"] += " "
os.environ["CFLAGS"] += ' '.join(compiler_flags)
base_path = sys.prefix

utils_dependence = ["cpp/src/utils.cpp"]
specialfunc_dependence = ["cpp/src/specialfunc.cpp"]
kerr_dependence = ["cpp/src/kerr.cpp", *specialfunc_dependence]
geo_dependence = ["cpp/src/geo.cpp", *specialfunc_dependence]
cf_dependence = ["cpp/src/cf.cpp", *utils_dependence]
swsh_dependence = ["cpp/src/swsh.cpp", *specialfunc_dependence]
monodromy_dependence = ["cpp/src/monodromy.cpp", *specialfunc_dependence]
hyperf_dependence = ["cpp/src/hypergeo_f.cpp", *specialfunc_dependence]
bessel_dependence = ["cpp/src/bessel.cpp", *hyperf_dependence, *specialfunc_dependence]
hyperu_dependence = ["cpp/src/hypergeo_u.cpp", *cf_dependence, *bessel_dependence, *specialfunc_dependence]
nusolver_dependence = ["cpp/src/nusolver.cpp", *monodromy_dependence, *cf_dependence, *swsh_dependence]
mst_dependence = ["cpp/src/mst.cpp", *nusolver_dependence, *hyperu_dependence, *hyperf_dependence]
gsn_dependence = ["cpp/src/gsn_asymp.cpp", *specialfunc_dependence]
radial_dependence = ["cpp/src/radialsolver.cpp", *mst_dependence, *gsn_dependence]
sourceint_dependence = ["cpp/src/sourceintegration.cpp", *geo_dependence, *swsh_dependence, *radial_dependence, "cpp/src/specialfunc.cpp"]
teuk_dependence = ["cpp/src/teukolsky.cpp", *sourceint_dependence, *radial_dependence]
hertz_dependence = ["cpp/src/hertz.cpp", *teuk_dependence]
metriccoeffs_dependence = ["cpp/src/metriccoeffs.cpp", *geo_dependence, *utils_dependence]
fluxes_dependence = ["cpp/src/fluxes.cpp", *teuk_dependence]
metric_dependence = ["cpp/src/metric.cpp", *hertz_dependence, *utils_dependence, *metriccoeffs_dependence, *kerr_dependence]
redshift_dependence = ["cpp/src/redshift.cpp", *metric_dependence]
unit_dependence = ["cpp/src/unit_test.cpp", *fluxes_dependence]

full_dependence = [*geo_dependence, *teuk_dependence, *radial_dependence, *hertz_dependence, *metriccoeffs_dependence, *fluxes_dependence, *redshift_dependence, *unit_dependence]

cpu_extension = dict(
    libraries=libraries,
    language='c++',
    include_dirs=["cpp/include", np.get_include()],
)

teuk_ext = Extension(
    "cybhpt_full", 
    sources=["cython/redshift_wrap.pyx", *set(full_dependence)], 
    **cpu_extension,
)

ext_modules = [teuk_ext]

setup(
    name="pybhpt",
    author="Zach Nasipak",
    author_email="znasipak@gmail.com",
    version = "0.0.1",
    description = "Black Hole Perturbation Theory and Self-Force Algorithms in Python",
    ext_modules = cythonize(ext_modules, language_level = "3"),
    packages=["pybhpt"],
    py_modules=["pybhpt.geo", "pybhpt.swsh", "pybhpt.radial", "pybhpt.teuk", "pybhpt.hertz", "pybhpt.metric", "pybhpt.flux"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Programming Language :: C++",
        "Programming Language :: Cython",
    ],
    cmdclass = {'build_ext': build_ext},
    zip_safe=False
)