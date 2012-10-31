#!/bin/bash

for f in $@; do
  echo pdfcrop $f $f
  pdfcrop $f $f
done

