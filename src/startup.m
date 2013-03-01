addpath(pwd());
addpath(fullfile(pwd, 'prob'));
addpath(fullfile(pwd, 'experiment'));
addpath(fullfile(pwd, 'plot'));
addpath(fullfile(pwd, 'util'));
addpath(fullfile(pwd, 'radonLikeFeatures'));

% DDP
addpath(genpath(fullfile(pwd, 'ddp')));

DATA_PREFIX = fullfile(pwd, '../data/kasthuri11/');

% Add local overrides here (i.e. cluster computation)
local_startup;

S1 = load_slices(fullfile(DATA_PREFIX, 'kasthuri11_mid.hd5'), 1, 1);
%A1 = load_ann_slices('../data/kasthuri11/kasthuri11_ann_mid.hd5', 1, 1);
%f1 = S1(1000:1199,1000:1199);
%f1 = S1(1200:1399,900:1099);
%f1 = S1(1400:1999,200:799);
f1 = S1(1600:1899,1:300);
%a1 = A1(1000:1199,1000:1199);
[PHI, misc_data] = feature_map(f1);

