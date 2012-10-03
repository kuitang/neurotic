%% Load from HD5
bock11 = load_slices('../data/bock11/bock11_begin.hd5', 1, 2);

%% Compute SIFT keypoints
I1 = bock11(:,:,1);
I2 = bock11(:,:,2);
[f1, d1] = vl_sift(bock11(:,:,1));
[f2, d2] = vl_sift(bock11(:,:,2));

%% Slow part
[matches, scores] = vl_ubcmatch(d1, d2);

%% Plots
[~, match_idx] = plot_sift(I1, I2, d1, d2, f1, f2, matches, scores);

Nmatch = length(scores);
df = f1(1:2, matches(1, match_idx)) - f2(1:2, matches(2, match_idx));
df = sqrt(sum(df .^ 2));

figure
hist(df);
title('Euclidean distances of matched keypoints');


addpath([docroot '/techdoc/creating_plots/examples'])
figure
h1 = vl_plotframe(f1(:,matches(1,match_idx)));
h2 = vl_plotframe(f2(:,matches(2,match_idx)));
set(h1, 'color', 'r', 'linewidth', 0.2);
set(h2, 'color', 'g', 'linewidth', 0.2);
title('Matched keypoints (red = I1, green = I2)');

% for i = 1 : N
%     m1 = matches(1, match_idx(i));
%     m2 = matches(2, match_idx(i));
%     [x, y] = dsxy2figxy([f1(1,m1) f1(2,m1)], ...
%                         [f2(1,m2) f2(2,m2)]);    
%     annotation('arrow', x, y);
% end

%%
figure
hist(scores)