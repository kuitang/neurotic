#!/bin/bash

. secret.m

URL="http://openconnecto.me/emca/$KAT11_TOKEN/objects/voxels/1"

for f in $(ls queries/); do
  prefix=$(echo $f | cut -d '_' -f 3)
  suffix=$(echo $f | cut -d '_' -f 4)

  outfile="kat11_ann_${prefix}_$suffix"

  echo /usr/bin/curl -o $outfile --data-binary @queries/$f $URL
done

