#!/usr/bin/bash
#module purge
#module load serverproxy

#rm -rf .eggs here build *.egg-info dist
# topdir=$PWD

#
# NOTE: Either source this script by 'source install.sh'
#       or execute it by './install.sh' and set the following paths manually.
# export PYTHONUSERBASE=${topdir}/here
# export PATH=${PYTHONUSERBASE}/bin:$PATH
# export PYTHONPATH=${PYTHONUSERBASE}/lib/python3.6/site-packages:${PYTHONPATH}

# install Python3 prerequisites
# pip3 install --ignore-installed --user pip
# pip install --user numpy==1.18.5
# pip install --user Cython
# pip install --user wheel

# build and install CMake
# git clone -b v3.23.2 https://github.com/Kitware/CMake.git build/cmake
# mkdir -p "${topdir}/build/cmake/build"
# cd "${topdir}/build/cmake/build"
# ../bootstrap --verbose --parallel=10 --prefix="${PYTHONUSERBASE}" -- "-DCMAKE_BUILD_TYPE:STRING=Release"
# make -j
# make install
# cd "${topdir}"

# build and install Swig - RHEL7.9 version doesn't work for openbabel
# git clone -b v3.0.12 https://github.com/swig/swig.git build/swig
# mkdir -p "${topdir}/build/swig/build"
# cd "${topdir}/build/swig"
# ./autogen.sh
# cd "${topdir}/build/swig/build"
# ../configure "--prefix=${PYTHONUSERBASE}"
# make -j
# make install
# cd "${topdir}"

# build and install openbabel with Python3 bindings
# git clone -b openbabel-3-1-1 https://github.com/openbabel/openbabel.git build/openbabel
# mkdir -p "${topdir}/build/openbabel/build"
# cd "${topdir}/build/openbabel/build"
# cmake "-DPYTHON_BINDINGS=ON" "-DPYTHON_EXECUTABLE=/usr/bin/python3" "-DRUN_SWIG=ON" "-DWITH_MAEPARSER=OFF" "-DCMAKE_INSTALL_PREFIX=${PYTHONUSERBASE}" ..
# sed -i 's|lib64|lib|' external/coordgen-master/coordgen/cmake_install.cmake
# sed -i 's|lib64|lib|' scripts/cmake_install.cmake
# make -j
# make install
# cd "${topdir}"

# install precomplex_generator
pip install .


# some tests
printf "NUMPY PATH: "
python3 -c 'import numpy as np; print(np.__file__)'

printf "PREC_GEN PATH: "
which prec_gen

