#!/bin/bash
# Download 100 slices from the beginning, middle, and end at full resolution.

# format is /res/xlow,xhigh/ylow,yhigh/zlow,zhigh/'
URLBASE='http://openconnecto.me/emca/bock11/hdf5/0/50000,52000/50000,52000'

/usr/bin/curl -o bock11_begin.hd5 $URLBASE/2917,3017/
/usr/bin/curl -o bock11_mid.hd5 $URLBASE/3500,3600/
/usr/bin/curl -o bock11_end.hd5 $URLBASE/4050,4150/

