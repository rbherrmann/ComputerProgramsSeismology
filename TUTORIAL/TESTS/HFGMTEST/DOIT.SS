#!/bin/sh

GDIR=${GREENDIR}/nCUSHF.REG
HS=0090
H=9.0
STK=5
DIP=80
RAKE=5
MW=3.0

while read DIST DELTA NSAMP T0 VRED D PROTO
do
rm -f ${PROTO}.Z
rm -f ${PROTO}.R
rm -f ${PROTO}.T
#for AZ in 10 
for AZ in 10 20 30 40 50 60 70 80 90 \
	100 110 120 130 140 150 160 170 180 \
	190 200 210 220 230 240 250 260 270 280 \
	290 300 310 320 330 340 350
do
gsac << EOF
mt to ZRT AZ ${AZ} STK ${STK} DIP ${DIP} RAKE ${RAKE} MW ${MW} FILE ${GDIR}/${HS}/${PROTO}
triangle width 0.5
w
q
EOF
saclhdr -NL -DIST -AZ -DEPMAX -DEPMIN T.Z >> ${PROTO}.Z
saclhdr -NL -DIST -AZ -DEPMAX -DEPMIN T.R >> ${PROTO}.R
saclhdr -NL -DIST -AZ -DEPMAX -DEPMIN T.T >> ${PROTO}.T
done
done < ${GDIR}/${HS}/W.CTL

#####
#   clean up
#####
rm -f T.Z T.R T.T

#####
#   now plot up the peak values
# in terms of hypocentral distance
#####
cat >> awkscr << EOF
{ print sqrt(\$1 * \$1 + ${H}*${H} ), \$3}
EOF
rm -f [ZRT].dat
for C in Z R T
do
cat *.${C} | awk -f awkscr >> ${C}.dat
done

#####
# plot all
#####
cat > acmdfil << EOF
'Z.dat' 1000 0.000 'SQ' 0.02
'R.dat'    1 0.000 'TR' 0.02
'T.dat' 1100 0.000 'CI' 0.02
EOF

genplt -A acmdfil -YLOG -XLOG -YMIN 0.0000001 -YMAX 0.001 -XMIN 1 -XMAX 300 -TX 'Hypocentral Distance (km)' -TY 'Vmax (m/s)'

calplt << EOF
NEWPEN
1
CENTER
4.5 6.2 0.14 'Vmax: Z(red) R(black) T(blue)' 0.0
PEND
EOF

cat GENPLT.PLT CALPLT.PLT > ALL.PLT

#####
#  plot up individual components
#####
cat > acmdfil << EOF
'Z.dat' 1000 0.000 'SQ' 0.02
EOF

genplt -A acmdfil -YLOG -XLOG -YMIN 0.0000001 -YMAX 0.001 -XMIN 1 -XMAX 300 -TX 'Hypocentral Distance (km)' -TY 'Vmax (m/s)'

calplt << EOF
NEWPEN
1
CENTER
4.5 6.2 0.14 'Vmax: Z(red) ' 0.0
PEND
EOF

cat GENPLT.PLT CALPLT.PLT > Z.PLT

#####
#  plot up individual components
#####
cat > acmdfil << EOF
'R.dat'    1 0.000 'TR' 0.02
EOF

genplt -A acmdfil -YLOG -XLOG -YMIN 0.0000001 -YMAX 0.001 -XMIN 1 -XMAX 300 -TX 'Hypocentral Distance (km)' -TY 'Vmax (m/s)'

calplt << EOF
NEWPEN
1
CENTER
4.5 6.2 0.14 'Vmax: R(black) ' 0.0
PEND
EOF

cat GENPLT.PLT CALPLT.PLT > R.PLT

#####
#  plot up individual components
#####
cat > acmdfil << EOF
'T.dat' 1100 0.000 'CI' 0.02
EOF

genplt -A acmdfil -YLOG -XLOG -YMIN 0.0000001 -YMAX 0.001 -XMIN 1 -XMAX 300 -TX 'Hypocentral Distance (km)' -TY 'Vmax (m/s)'

calplt << EOF
NEWPEN
1
CENTER
4.5 6.2 0.14 'Vmax: T(blue) ' 0.0
PEND
EOF

cat GENPLT.PLT CALPLT.PLT > T.PLT

