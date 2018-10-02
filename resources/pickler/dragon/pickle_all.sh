#!/bin/bash

#### 
#### By Elim Thompson (09/25/2018)
####
#### This is a bash script that runs all
#### pickling python scripts.
####
##############################################

#### step 1. neutrinos
for dtype in numu nue nutau; do
    echo "======================================="
    echo "pickling ${dtype} ..."
    python genie_nu.py --flavor ${dtype}
done

#### step 2. muon
echo "======================================="
echo "pickling muon ..."
python muon.py

#### step 3. noise
#echo "======================================="
#echo "pickling noise ..."
#python noise.py

#### step 4. data
echo "======================================="
echo "pickling data ..."
python data.py

#### bye!
echo "======================================="
echo "bye !"
