addpath(fullfile(pwd, '.'));
addpath(fullfile(pwd, 'prob'));
addpath(fullfile(pwd, 'experiment'));
addpath(fullfile(pwd, 'plot'));
addpath(fullfile(pwd, 'util'));

DATA_PREFIX = fullfile(pwd, '../data/kasthuri11/');

% Add local overrides here (i.e. cluster computation)
local_startup;

S1 = load_slices(fullfile(DATA_PREFIX, 'kasthuri11_mid.hd5'), 1, 1);
%A1 = load_ann_slices('../data/kasthuri11/kasthuri11_ann_mid.hd5', 1, 1);
f1 = S1(1000:1199,1000:1199);
%a1 = A1(1000:1199,1000:1199);
PHI1 = feature_map(f1);

