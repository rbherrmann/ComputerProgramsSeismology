#!/bin/sh
export PATH=:/home/rbh/PROGRAMS.310t/SPHERICITY/MINEOS/bin:$PATH
#
# Example of the synthetic seismogram's generation script.
#
# Usage: RUN_MINEOS.sh model_name
# 
# Available model names for DEMO version:
# prem_noocean  prem_ocean CPacific NRussia tak135-f
#
#=========================================================
#
# function creat_orig
#
creat_origin() {
time=`awk '{if(NR == 1){yd=$2-1970;vis=int((yd+1)/4); s=(365*yd+vis+$3-1)*86400.0+($4*60+$5)*60+$6; printf("%17.5f", s);}}' < $1`
awk 'BEGIN{t='"$time"';}{if(NR == 1) \
printf("%9.4f %9.4f %9.4f %17.5f %8d %8d %8d %4d %4d %4d %8d %8d %-7.7s %9.4f %-1.1s %7.2f %8d %7.2f %8d %7.2f %8d %-15.15s %-15.14s %8d %-17.17s\n", \
$7,$8,$9,t,1,-1,$2*1000+$3,-1,-1,-1,-1,-1,"-",-999.0000,"-",-999.0,-1,-999.0,-1, \
-999.00,-1,"-","PDE & Hvd CMT",-1,-1); \
}' $1 > $2.origin ;
}
#
#   Main procedure
#
if test "$#" != 1; then
echo " Usage: RUN_MINEOS.sh model_name"
exit
fi
model=$1                # setup 1-D model name
# check model name
flg=0
for f in \
CPacific NRussia tak135-f prem_noocean prem_noocean_na prem_ocean prem_noocean_1ln.txt model.dk nnak135.coef noqtak135-f
do
if test "$f" = $1; then
flg=1
fi
done
if test "$flg" = 0; then
echo "Model name $1 is wrong, allowed names are:"
echo "CPacific NRussia tak135-f prem_noocean prem_noocean_na prem_ocean prem_noocean_1ln.txt model.dk"
exit
fi
#=========================================================
# 1. run minos_bran program for fundamental S  mode,
# where,  n=0, 0 < f <  0.2 Hz,
#
echo "Step 1:  minos_bran runs for S modes ....................."
echo "============== Program minos_bran =================="
null=`ls * | grep '_S$'`
if test "X$null" != X; then 
rm -f $null
fi
time minos_bran << EOF
../models/${model}.txt
${model}_S
e${model}_S
1.0e-10 10
3
0 2000 0.0 200.0 0 3
EOF
#=========================================================
# 1a. run minos_bran program for fundamental R  mode,
# where,  n=0, 0 < f <  0.2 Hz,
#
#echo "Step 1a:  minos_bran runs for R modes ....................."
#echo "============== Program minos_bran =================="
#null=`ls * | grep '_R$'`
#if test "X$null" != X; then 
#rm -f $null
#fi
#time minos_bran << EOF
#../models/${model}.txt
#${model}_R
#e${model}_R
#1.0e-10 1
#1
#2 8000 0.0 200.0 0 2000
#EOF
#=========================================================
# 2. run minos_bran program for fundamental T  mode,
# where,  n=0, 0 < f <  0.2 Hz,
#
echo "Step 2: minos_bran runs for T modes ....................."
echo "============== Program minos_bran =================="
null=`ls * | grep '_T$'`
if test "X$null" != X; then 
rm -f $null
fi
time minos_bran << EOF
../models/${model}.txt
${model}_T
e${model}_T
1.0e-10 10
2
0 2000 0.0 200.0 0 3
EOF
# 3. Convert minos_bran results to .eigen relation (S mode)
#
echo "Step 3: eigen for S ....................................."
if test -f test_S.eigen; then
 rm -rf test_S.*
fi
time eigcon << EOF
3
../models/${model}.txt
1000
${model}_S
e${model}_S
test_S
EOF
#============================================================
# 3a. Convert minos_bran results to .eigen relation (R mode)
#
#echo "Step 3a: eigen for R ....................................."
#if test -f test_R.eigen; then
# rm -rf test_R.*
#fi
#time eigcon << EOF
#3
#../models/${model}.txt
#1000
#${model}_R
#e${model}_R
#test_R
#EOF
#============================================================
# 4. Convert minos_bran results to .eigen relation (T mode)
#echo "Step 4: eigen for T ....................................."
#if test -f test_T.eigen; then
 #rm -rf test_T.*
#fi
#time eigcon << EOF
#2
#../models/${model}.txt
#1000
#${model}_T
#e${model}_T
#test_T
#EOF
#RBH_#=========================================================
#RBH_# 5. Evaluate green functions for given sitechan relation
#RBH_echo "Step 5: green functions evaluation ........................."
#RBH_if test -f green.wfdisc; then
#RBH_ rm -rf green.*
#RBH_fi
#RBH_time green << EOF
#RBH_short
#RBH_db_list
#RBH_china_cmt_event
#RBH_1 260
#RBH_8000
#RBH_green
#RBH_EOF
#RBH_cp -p short.site green.site
#RBH_cp -p short.sitechan green.sitechan
#RBH_# create origin relation for data base green
#RBH_../scripts/creat_origin china_cmt_event green
#RBH_#============================================================
#RBH_# 6. Synthetic data construction
#RBH_echo "Step 6: synthetic seismogram construction .................."
#RBH_if test -f Syndat.wfdisc; then
#RBH_ rm -rf Syndat.*
#RBH_fi
#RBH_time syndat << EOF
#RBH_china_cmt_event
#RBH_0
#RBH_green
#RBH_Syndat
#RBH_0
#RBH_EOF
#RBH_cp -p short.site Syndat.site
#RBH_cp -p short.sitechan Syndat.sitechan
#RBH_creat_origin china_cmt_event Syndat
#RBH_
#RBH_#####
#RBH_#	step 5 - make SAC files
#RBH_#####
#RBH_cucss2sac Syndat Syndat_SAC
#RBH_
#RBH_
