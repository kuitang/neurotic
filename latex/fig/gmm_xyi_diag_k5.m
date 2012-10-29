% Revision 2cc0c7b
load('gmm_xyi_diag_k5_data.mat')
idxs = randsample(40000,10000);
imshow(f1_clahe);
hold on;
gscatter(PHI(idxs,1), PHI(idxs,2), B.z(idxs));
title('GMM (x,y,i), diag, K=5');
