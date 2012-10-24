#!/bin/bash
# Download 100 slices from the beginning, middle, and end at full resolution.

# format is /res/xlow,xhigh/ylow,yhigh/zlow,zhigh/'
URLBASE='http://openconnecto.me/emca/kasthuri11/hdf5/0/12000,14000/12000,14000/'

/usr/bin/curl -o kasthuri11_begin.hd5 $URLBASE/400,500/
/usr/bin/curl -o kasthuri11_mid.hd5 $URLBASE/900,1000/
/usr/bin/curl -o kasthuri11_end.hd5 $URLBASE/1750,1850

