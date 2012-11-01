function [ slices ] = load_slices( hdffile, first, last)
% load_slices loads single-precision [0, 1] images from openconnecto.me
%   slices = load_slices(hdffile, idxs) gives a 3-dimensional array of the
%   slices given in idxs from hdffile, format given by openconnecto.me
%   volume cut service.

    f = h5read(hdffile, '/cube', [1 1 first], [Inf Inf last]);
    slices = single(f) / 255;

end

