#!/bin/sh


#####
#       define the color palete for elevation
#       MIN_ELV R G B MAX_ELEV R G B
#       (from Chuck Ammon)
#####
cat > map.cpt << EOF
-1000   100     200     255     -500    100     200     255
-500    150     225     255     0       150     225     255
0       100     150     100     30      100     150     100
30      125     175     125     60      125     175     125
60      150     200     150     122     150     200     150
122     175     225     175     183     175     225     175
183     200     255     200     244     200     255     200
244     212     255     212     305     212     255     212
305     255     255     225     457     255     255     225
457     255     225     175     610     255     225     175
610     255     225     125     702     255     225     125
702     255     175     75      914     255     175     75
914     200     150     50      1219    200     150     50
1219    175     125     50      1450    175     125     50
1450    150     100     50      1700    150     100     50
1700    150     125     100     1981    150     125     100
1981    125     125     125     2134    125     125     125
2134    150     150     150     2438    150     150     150
2438    175     175     175     2743    175     175     175
2743    200     200     200     3048    200     200     200
3048    233     233     233     9000    233     233     233
B       100     200     255
F       100     200     255
EOF

#####
#    set GMT defaults
#####
gmtset BASEMAP_TYPE FANCY DEGREE_FORMAT 1 MEASURE_UNIT cm 

#####
#    Define the default name of the PostScript output
#####
FNAME="map.eps"


#####
#    Define raster for topography:  
#####
GRDRAS="1"

#####
#    Define map bounds: MINLON/MAXLON/MINLAT/MAXLAT 
#####
LATLON="-4.000000/0.000000/35.000000/39.000000"

#####
#    Define Mercalli projection: Center_lon/Center_Lat/Plot_Width   
#####
PROJ="M-2.000000/37.000000/15c"
rm -f ${FNAME}

#####
#    Define Coastline resolution: one of clihf 
#####
RESCOAST="i"

#####
#    Define Ticmark interval 
#####
TICS="a1.000000g0/a1.000000g0WSne"

#####
#    Define epicenter symbol and size  
#####
EPISYM="A0.50c"
EPISYM="c0.37c"
EPICOLOR="255/0/0"

#####
#    Define station symbol and size  
#####
STASYM="c0.37c"
STACOLOR="255/0/0"

#####
#    Define boundaries for pscoast 
#    1 = National boundaries 
#    2 = State boundaries within the Americas 
#    3 = Marine boundaries 
#    a = All boundaries (1-3) pscoast 
#####
BDRYS="a"

#####
#    Define resolution for database raster assuming 5m for GRDRAS=1  
#####
RESIN="5m"

#####
#    Define resolution for output of grdsample 
#####
RESOUT="2m"

grdraster ${GRDRAS} -I${RESIN} -R${LATLON} -Gcmap.grd -V
grdsample cmap.grd -Gmap.grd -I${RESOUT} -R${LATLON} -V
grdgradient map.grd -A135 -GMA.grd -Nt -V 
grdimage -P map.grd -X2.5i -Y1.5i -J${PROJ} -R${LATLON} -Cmap.cpt -K -IMA.grd  -V > ${FNAME} 

#####
#    Plot coast 
#####
pscoast -P -W  -J${PROJ} -R${LATLON} -B${TICS} -O  -K  -N${BDRYS} -D${RESCOAST} -A2500  -V  >> ${FNAME}
project -C-2.582200/37.189701 -E-2.005700/37.039398 -G10 -Q |
	psxy -P -R${LATLON} -J${PROJ} -W2.0 -O     -K >> ${FNAME}
project -C-2.582200/37.189701 -E-1.988000/37.583801 -G10 -Q |
	psxy -P -R${LATLON} -J${PROJ} -W2.0 -O     -K >> ${FNAME}
project -C-2.005700/37.039398 -E-1.988000/37.583801 -G10 -Q |
	psxy -P -R${LATLON} -J${PROJ} -W2.0 -O     -K >> ${FNAME}

#####
#    Plot epicenters 
#####
psxy -P -J${PROJ} -R${LATLON} -O -: -S${EPISYM} -W0.8 -G${EPICOLOR} -V -K   >> ${FNAME} << EOF
37.189701 -2.582200
37.189701 -2.582200
37.039398 -2.005700
EOF

#####
#    Plot stations 
#####
psxy -P -J${PROJ} -R${LATLON} -O -: -S${STASYM} -W0.8 -G${STACOLOR} -V   >> ${FNAME} << EOF
37.039398 -2.005700
37.583801 -1.988000
37.583801 -1.988000
EOF


######
#     Cleanup
######
rm -f map.grd cmap.grd MA.grd map.cpt
