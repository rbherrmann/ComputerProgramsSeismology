#!/bin/csh
#
# Example of the BIG_ENDIAN <===> LOW_ENDIAN conversion.
#
# Usage: RUN_ENDIAN.csh db_name
#
#=========================================================
alias cp cp
alias rm rm
alias mv mv
#=========================================================
if( $#argv != 1) then
echo " Usage: RUN_ENDIAN.csh db_name"
exit
endif
if ( ! -f $1.wfdisc) then
echo " .wfdisc relation $1.wfdisc does not exist"
exit
endif
if ( ! -d $1.wfdisc.dat) then
echo " Directory $1.wfdisc.dat does not exist"
exit
endif
#
set big=`egrep -e"t4" $1.wfdisc | wc -l`
set low=`egrep -e"f4" $1.wfdisc | wc -l`
if($big != 0 && $low != 0) then
echo ".wfdisc inclides BIG and LOW orders"
exit
endif
cd $1.wfdisc.dat
if($big != 0) then
echo " BIG_ENDIAN ==> LOW_ENDIAN conversion"
sed -e"s/t4/f4/" ../$1.wfdisc > tmp.$$
mv tmp.$$ ../$1.wfdisc
endi 4 [wg].*
else
echo " LOW_ENDIAN ==> BIG_ENDIAN conversion"
sed -e"s/f4/t4/" ../$1.wfdisc > tmp.$$
mv tmp.$$ ../$1.wfdisc
endi 4 [wg].*
endif
cd ..
exit
