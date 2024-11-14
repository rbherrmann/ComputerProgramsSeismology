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
gmtset BASEMAP_TYPE FANCY DEGREE_FORMAT 5 MEASURE_UNIT cm 

#####
#    Define the default name of the PostScript output
#####
FNAME="map.eps"


#####
#    Define map bounds: MINLON/MAXLON/MINLAT/MAXLAT 
#####
LATLON="-80.811996/-60.266800/-56.453602/-33.860100"

#####
#    Define raster for topography:  
#####
GRDRAS="1"

#####
#    Define Mercalli projection: Center_lon/Center_Lat/Plot_Width   
#####
PROJ="M-71.039398/-45.156853/15c"
rm -f ${FNAME}

#####
#    Define Coastline resolution: one of fhilc 
#    (f)ull, (h)igh, (i)ntermediate, (l)ow, and (c)rude)
#  . The resolution drops off by 80% between data sets.
#####
RESCOAST="l"

#####
#    Define Ticmark interval 
#####
TICS="a2.000000g0/a5.000000g0WSne"

#####
#    Define epicenter symbol and size  
#####
EPISYM="A0.20"
EPICOLOR="255/255/0"

#####
#    Define station symbol and size  
#####
STASYM="c0.15"
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
#    Define resolution for output of grdsample 
#####
RESOUT="10m"

pscoast -P -W  -J${PROJ} -R${LATLON} -B${TICS} -G200  -K  -N${BDRYS} -D${RESCOAST} -A2500  -V   > ${FNAME}

#####
#    MAP GREAT CIRCLE RAY PATHS USING PROJECT
#####
psxy -P -J${PROJ} -R${LATLON} -W2.0 -O -: -V -K >> ${FNAME} << EOF
-36.360100 -73.811996 
-45.572990 -72.081390 
EOF
psxy -P -J${PROJ} -R${LATLON} -W2.0 -O -: -V -K >> ${FNAME} << EOF
-36.360100 -73.811996 
-53.953602 -68.266800 
EOF

#####
#    PLOT EPICENTER LOCATIONS
#####
psxy -P -J${PROJ} -R${LATLON} -O -: -S${EPISYM} -W0.8 -G${EPICOLOR} -V -K   >> ${FNAME} << EOF
-36.360100 -73.811996
-36.360100 -73.811996
EOF


#####
#    PLOT STATION LOCATIONS
#####
psxy -P -J${PROJ} -R${LATLON} -O -: -S${STASYM} -W0.8 -G${STACOLOR} -V -K >> ${FNAME} << EOF
-45.572990 -72.081390
-53.953602 -68.266800
EOF

pstext -P -J${PROJ} -R${LATLON} -O -:  -G0/0/0 -V -Dj0.1i/0.1i  >> ${FNAME} << EOF
-45.572990 -72.081390 14 0 1 LB COYC
-53.953602 -68.266800 14 0 1 LB DSPA
EOF


######
#     Cleanup
######
rm -f map.grd cmap.grd MA.grd map.cpt
