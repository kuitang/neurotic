function [ gmm ] = make_gmm_prior( X, K, background_like )
    [N, D] = size(X);
    nX = max(X(:,1));
    nY = max(X(:,2));
    
    % Want a STRONG pri for a HIGH variance (LOW precision).
    % So TINY gamma.
    % Got the settings from this Wikipedia page..    
    gmm = struct('K', K, ...
                 'N', N, ...
                 'nX', nX, ...
                 'nY', nY, ...
                 'prior_mean', [100 100 0.75], ...
                 'prior_cov', cov(X), ...
                 'prior_dof', 4, ...
                 'prior_scale', 0.1, ... % Contribute more precision??? UNDERSTAND THIS PARAMETER!
                 's_z', randsample(K, N, true), ...
                 'prior_mix', []);        
    
    % Weaker Dirichlet prior (\sum \alpha_k = 1)
    gmm.prior_mix = zeros(K, 1);
    gmm.prior_mix(1) = 0.4;
    gmm.prior_mix(2:end) = 0.6 / (K - 1);        
    
    % All points with intensity < 0.3 are considered candidates for
    % determined by a very sophisticated method called "taking thresholds
    % and looking." For now, fit a method-of-moments Beta (notebook pp 23)    
        
    gmm.background_like = background_like;        
    
    gmm.h_diagnostic = figure;
    set(gmm.h_diagnostic, 'Units', 'normalized', 'Position', ...
        [0 0.7 0.9 0.4]);
    gmm.plot_diagnostic = @plot_intensity_hists;
    
    % Initialize with k-means
    % Only valid if each feature is Gaussian!
    % [gmm.z, gmm.mu_mean] = kmeans(X, K);
        
end
