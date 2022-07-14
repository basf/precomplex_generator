#!/usr/bin/bash

#rm -rf .eggs install build *.egg-info dist

export PYTHONUSERBASE=$PWD/install

pip3 install --ignore-installed --user pip
export PATH=$PYTHONUSERBASE/bin:$PATH
which pip

pip install --user numpy==1.18.5
pip install --user Cython
pip install --user wheel
pip install --user .
python3 setup.py install --user 

export PYTHONPATH=$PYTHONUSERBASE/lib/python3.6/site-packages:$PYTHONPATH

#printf "NUMPY PATH: "
#python3 -c 'import numpy as np; print(np.__file__)'

#printf "PREC_GEN PATH: "
#which prec_gen
