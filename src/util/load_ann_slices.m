function [ ann_slices ] = load_ann_slices( hdffile, first, last)
% load_slices loads uint32 identifiers from openconnecto.me

    ann_slices = h5read(hdffile, '/cube', [1 1 first], [Inf Inf last]);

end

