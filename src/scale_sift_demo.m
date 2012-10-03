%% Load from HD5
bock11 = load_slices('../data/bock11/bock11_begin.hd5', 1, 2);

%% Compute SIFT keypoints
I1 = bock11(:,:,1);
I2 = bock11(:,:,2);
[f1, d1] = vl_sift(bock11(:,:,1));
[f2, d2] = vl_sift(bock11(:,:,2));

%% cluster
sc1 = f1(3, 1:(0.9*end));
sc2 = f2(3, 1:(0.9*end));

all_sc = [sc1 sc2];
all_idxs = kmeans(all_sc, 5);

idxs1 = all_idxs(1:length(sc1));
idxs2 = all_idxs((length(sc1) + 1):end);

%% match
matches = [];
scores = [];
for k = 1 : 5
    dd1 = d1(:, idxs1(idxs1 == k));
    dd2 = d2(:, idxs2(idxs2 == k));
    [m, s] = vl_ubcmatch(dd1, dd2);
    matches = [matches m];
    scores = [scores s];
end

%% ?
matches
scores

%% plot
plot_sift(I1, I2, matches, scores);
