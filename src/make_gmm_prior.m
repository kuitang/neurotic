function [ gmm ] = make_gmm_prior( X )
    [N, D] = size(X);
    nX = max(X(:,1));
    nY = max(X(:,2));
    
    % Want a STRONG pri for a HIGH variance (LOW precision).
    % So TINY gamma.
    % Got the settings from this Wikipedia page..    
    gmm = struct('K', 6, ... % To begin, only the background class.
                 'Kmax', 500, ....
                 'N', N, ...
                 'nX', nX, ...
                 'nY', nY, ...
                 's_z', zeros(N, 1), ...
                 'prior_mean', [100 100 0.75], ...
                 'prior_cov', cov(X), ...
                 'prior_dof', 5, ...
                 'prior_scale', 0.5, ... % Number of prior observations
                 'prior_conc', 1);
    
    % Tighter intensities
    %gmm.prior_cov(3,3) = 0.01 * gmm.prior_cov(3,3);             
    
    % Precompute predictive likelihoods for zero-point clusters    
    
    v_intensity = gmm.prior_cov(3,3);
    gmm.prior_cov(3,:) = zeros(3, 1);
    gmm.prior_cov(:,3) = zeros(1, 3);
    gmm.prior_cov(3,3) = 0.1 * v_intensity;    
    
    % Weaker Dirichlet prior (\sum \alpha_k = 1)
%     gmm.prior_mix = zeros(K, 1);
%     gmm.prior_mix(1) = 0.4;
%     gmm.prior_mix(2:end) = 0.6 / (K - 1);
%     
%     gmm.prior_mix = 100 * gmm.prior_mix;    
       
    gmm.background_pdf = make_sigmoid_pdf(0.5, 30);        
    
    gmm.s_z = kmeans(X, gmm.K - 1);
    % shift up
    gmm.s_z = gmm.s_z + 1;
    
    bg = rand < gmm.background_pdf(X(:,3));
    % assign to background
    gmm.s_z(bg) = 1;    
    
    gmm.h_diagnostic = figure;
    set(gmm.h_diagnostic, 'Units', 'normalized', 'Position', ...
        [0 0.7 1 0.4]);
    %gmm.plot_diagnostic = @plot_intensity_hists;
    
    % Initialize with k-means, but don't set background
    % Only valid if each feature is Gaussian!    
    %[gmm.s_z, ~] = kmeans(X, K - 1);
    %gmm.s_z = gmm.s_z + 1;
        
end
