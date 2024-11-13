#!/bin/csh
#
# Example of the synthetic seismogram's generation script.
#
# Usage: RUN_MINEOS.csh model_name
# 
# Available model names for DEMO version:
# prem_noocean  prem_ocean CPacific NRussia
#
#=========================================================
if( $#argv != 1) then
echo " Usage: RUN_MINEOS.csh model_name"
exit
endif
alias cp cp
alias rm rm
alias mv mv
set model=$1                # setup 1-D model name
# check model name
set flg=0
foreach f (CPacific NRussia prem_noocean prem_noocean_na prem_ocean prem_noocean_1ln.txt)
if($f == $1) set flg=1
end
if($flg == 0) then
echo "Model name $1 is wrong, allowed names are:"
echo "CPacific NRussia prem_noocean prem_noocean_na prem_ocean prem_noocean_1ln.txt"
exit
endif
#=========================================================
# 1. run minos_bran program for fundamental S  mode,
# where,  n=0, 0 < f <  0.2 Hz,
#
echo "Step 1:  minos_bran runs for S modes ....................."
echo "============== Program minos_bran =================="
if( -f ${model}_S) rm -f ${model}_S
if( -f e${model}_S) rm -f e${model}_S
time minos_bran << EOF
../models/$model.txt
${model}_S
e${model}_S
1.0e-10 1
3
2 8000 0.0 200.0 0 0
EOF
#=========================================================
# 2. run minos_bran program for fundamental T  mode,
# where,  n=0, 0 < f <  0.2 Hz,
#
echo "Step 2: minos_bran runs for T modes ....................."
echo "============== Program minos_bran =================="
if( -f ${model}_T) rm -f ${model}_T
if( -f e${model}_T) rm -f e${model}_T
time minos_bran << EOF
../models/$model.txt
${model}_T
e${model}_T
1.0e-10 1
2
2 8000 0.0 200.0 0 0
EOF
#============================================================
# 3. Convert minos_bran results to .eigen relation (S mode)
#
echo "Step 3: eigen for S ....................................."
if( -f test_S.eigen) rm -rf test_S.*
time eigcon << EOF
3
../models/$model.txt
1000
${model}_S
e${model}_S
test_S
EOF
#============================================================
# 4. Convert minos_bran results to .eigen relation (T mode)
echo "Step 4: eigen for T ....................................."
if( -f test_T.eigen) rm -rf test_T.*
time eigcon << EOF
2
../models/$model.txt
1000
${model}_T
e${model}_T
test_T
EOF
#=========================================================
# 5. Evaluate green functions for given sitechan relation
echo "Step 5: green functions evaluation ........................."
if( -f green.wfdisc) rm -rf green.*
time green << EOF
short
db_list
china_cmt_event
10 260
8000
green
EOF
cp -p short.site green.site
cp -p short.sitechan green.sitechan
# create origin relation for data base green
creat_origin china_cmt_event green
#============================================================
# 6. Synthetic data construction
echo "Step 6: synthetic seismogram construction .................."
if( -f Syndat.wfdisc) rm -rf Syndat.*
time syndat << EOF
china_cmt_event
0
green
Syndat
0
EOF
# convert Green functions to BIG_ENDIAN binary numerical storage
# this is ANTELOPE requirement
# $bin/wfendi 0 0 Syndat Syndat
# copy .site, .sitechan relations into database Syndat
cp -p short.site Syndat.site
cp -p short.sitechan Syndat.sitechan
creat_origin china_cmt_event Syndat
exit
