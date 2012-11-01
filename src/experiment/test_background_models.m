% Background models -- assume intensity is dimension 3
triangle   = @(gmm_, X_) 2 / gmm_.N * (1 - X_(:,3));
beta_tight = @(gmm_, X_) 1 / gmm_.N * betapdf(X_(:,3), 1, 3);
beta_loose = @(gmm_, X_)  1 / gmm_.N * betapdf(X_(:,3), 1, 1.5);
ss_sigmoid = @(gmm_, X_) (1 / gmm_.N) * (1 / 0.21265926) * (1 - (1 ./ (1 + exp( -10 * (X_(:,3) - 0.2) ) ) ) );
ss_sigmoid_filter = @(gmm_, X_) (1 / gmm_.N) * (1 / 0.21265926) * (1 - (1 ./ (1 + exp( -10 * (X_(:,3) - 0.2) ) ) ) ) .* (background_filter(gmm_, X_, 1));



%% Run some samples
figure;
triangle_samps   = gmm_gibbs(PHI1, f1, make_gmm_prior(PHI1, 6, triangle), 500, 200, 10);
title('GMM triangular background, k = 6');
figure;
beta_tight_samps = gmm_gibbs(PHI1, f1, make_gmm_prior(PHI1, 6, beta_tight), 500, 200, 10);
title('GMM beta(1,3) background, k = 6');
figure;
beta_loose_samps = gmm_gibbs(PHI1, f1, make_gmm_prior(PHI1, 6, beta_tight), 500, 200, 10);
title('GMM beta(1,1.5) background, k = 6');
figure;
ss_sigmoid_samps = gmm_gibbs(PHI1, f1, make_gmm_prior(PHI1, 6, ss_sigmoid), 500, 200, 10);
title('GMM background sigmoid centered at 0.2, compressed 10x, k =6');

%% The shebang
figure;
ss_sigmoid_filter_samps = gmm_gibbs(PHI1, f1, make_gmm_prior(PHI1, 6, ss_sigmoid_filter), 500, 200, 10);
title('GMM background sigmoid center at 0.2, compressed 10x, 3x3 filter, k =6');

%% Collate samples
% triangle_mean   = struct_mean(triangle_samps);
% beta_tight_mean = struct_mean(beta_tight_samps);
% beta_loose_mean = struct_mean(beta_loose_samps);

%% Points of disagreement on background (take just the last sample...)
% Really need to get a good collater for MCMC...

% Compute the intersection of background points, and plot the deviation of
% each.
triangle_bg   = triangle_samps(end).s_z == 1;
beta_tight_bg = beta_tight_samps(end).s_z == 1;
beta_loose_bg = beta_loose_samps(end).s_z == 1;
ss_sigmoid_bg = ss_sigmoid_samps(end).s_z == 1;

bg_intersect = triangle_bg & beta_tight_bg & beta_loose_bg;
fprintf(1, 'Proportion background: triangle = %g beta tight = %g beta loose = %g intersection = %g\n', ...
        mean(triangle_bg), mean(beta_tight_bg), mean(beta_loose_bg), mean(bg_intersect));
    
triangle_d   = triangle_bg ~= bg_intersect;
beta_tight_d = beta_tight_bg ~= bg_intersect;
beta_loose_d = beta_loose_bg ~= bg_intersect;
fprintf(1, 'Disagreement with intersection: triangle = %g beta_tight = %g beta_loose = %g\n', ...
        mean(triangle_d), mean(beta_tight_d), mean(beta_loose_d));

fprintf(1, 'Disagreement of beta with each other: %g\n', mean(beta_tight_d ~= beta_loose_d));

%% Plot colors
c = eye(3);
[N, ~] = size(PHI1)
[Nx Ny] = size(f1);
dots = zeros(Nx, Ny, 3);
for n = 1:N
    ncolors = 0;
    x = PHI1(n,1);
    y = PHI1(n,2);
    if triangle_bg(n)
        dots(x,y,:) = squeeze(dots(x,y,:)) + c(1,:)';
        ncolors = ncolors + 1;
    end
    
    if beta_tight_bg(n)
        dots(x,y,:) = squeeze(dots(x,y,:)) + c(2,:)';
        ncolors = ncolors + 1;
    end
    
    if beta_loose_bg(n)
        dots(x,y,:) = squeeze(dots(x,y,:)) + c(3,:)';
        ncolors = ncolors + 1;
    end
    
    %assert(ncolors < 3);
    
    if ncolors > 0
        dots(x,y,:) = dots(x,y,:) / ncolors;            
    else
        dots(x,y,:) = [1 1 1]';
    end
    
end

figure
image(dots)
title('background distribution. R = triangle, B = beta(1,3), C = beta(1,1.5)')
