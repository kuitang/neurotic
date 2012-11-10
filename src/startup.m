addpath('.');
addpath('./prob');
addpath('./experiment');
addpath('./plot');
addpath('./util');

S1 = load_slices('../data/kasthuri11/kasthuri11_begin.hd5', 1, 1);
f1 = S1(1000:1199,1000:1199);
PHI1 = feature_map(f1);