%% Background models
triangle   = @(gmm_, X_) 2 / gmm_.N * (1 - X_(:,3));
ss_sigmoid = @(gmm_, X_) (1 / gmm_.N) * (1 / 0.21265926) * (1 - (1 ./ (1 + exp( -10 * (X_(:,3) - 0.2) ) ) ) );

%% Prior
mm = make_spatial_mm_prior(PHI1, 6);
mm.prior_beta_params = [ 0.1 0.1 0.0005 ];
mm.background_like = ss_sigmoid;
mm.beta_posterior_samps = cell(6, 1);
mm.feature_likes = { @beta_intensity_like };


%% Samples!
figure;
beta_samps = run_mm(@imm_iter, PHI1, f1, mm, 400, 100, 10);
title('IMM beta [0.1 0.1 0.005] 1-avoiding, k =6 ');