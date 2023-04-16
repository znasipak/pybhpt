import os
from Cython.Build import cythonize
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
from subprocess import run

compiler_flags = ["-std=c++11", "-Xclang", "-fopenmp", "-O2"]

output = run(["uname", "-m"], capture_output=True)
if output.stdout.decode("utf-8") == "x86_64":
    compiler_flags.append("-march=x86-64")

# I have to add this flag directly to the compiler because there
# does not seem to be any way to pass it to the compiler when 
# constructing the library
os.environ["CC"] += " "
os.environ["CC"] += ' '.join(compiler_flags)

utils_dependence = ["cpp/src/utils.cpp"]
specialfunc_dependence = ["cpp/src/specialfunc.cpp"]
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

full_dependence = [*geo_dependence, *teuk_dependence, *radial_dependence, *hertz_dependence, *metriccoeffs_dependence]

# create library of all the c++ files that get linked at the end in setup so that we are
# not duplicating c++ file compilation across the different extensions
lib_extension = dict(
    sources = [*set(full_dependence)],
    libraries=["gsl", "gomp", "boost_filesystem"],
    language='c++',
    include_dirs = ["cpp/include"],
    macros = None
)
sfcpp = ['sfcpp', lib_extension]

# libsfcpp_extension = Extension(
#     "libsfcpp",
#     sources = [*set(full_dependence)],
#     libraries=["gsl", "gomp", "boost_filesystem"],
#     extra_compile_args=["-std=c++11", "-Xclang", "-fopenmp"],
#     language='c++',
#     include_dirs = ["include"],
# )

cpu_extension = dict(
    libraries=["gsl", "gomp", "boost_filesystem"],
    language='c++',
    # extra_compile_args=["-Xclang", "-fopenmp", "-O2"], 
    # extra_compile_args=["-Xclang", "-fopenmp", "-O0", "-ggdb"], # use to speed up compilation during development
    include_dirs=["cpp/include", np.get_include()],
)

teuk_ext = Extension(
    "teukmodecy", 
    sources=["cython/teukolsky_wrap.pyx"], 
    **cpu_extension,
)

# ext_modules = [teuk_ext, geo_ext, radial_ext]
ext_modules = [teuk_ext]

setup(
    name="pybhpt",
    author="Zach Nasipak",
    author_email="znasipak@gmail.com",
    version = "0.0.1",
    description = "Self-Force Algorithms in Python",
    ext_modules = cythonize(ext_modules, language_level = "3"),
    packages=["pybhpt", "pybhpt.geo", "pybhpt.radial", "pybhpt.teuk", "pybhpt.hertz"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Programming Language :: C++",
        "Programming Language :: Cython",
    ],
    libraries = [sfcpp],
    cmdclass = {'build_ext': build_ext},
    zip_safe=False
)