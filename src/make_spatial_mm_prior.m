function [ gmm ] = make_spatial_mm_prior( X, K )
% Make a prior only for the spatial components.
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
                 'prior_mean', [100 100], ...
                 'prior_cov', cov(X(:,1:2)), ...
                 'prior_dof', 4, ...
                 'prior_scale', 1, ...
                 's_z', randsample(K, N, true), ...
                 'prior_mix', []);        
    
    % Weaker Dirichlet prior (\sum \alpha_k = 1)
    gmm.prior_mix = zeros(K, 1);
    gmm.prior_mix(1) = 0.4;
    gmm.prior_mix(2:end) = 0.6 / (K - 1);            
    
    % Initialize with k-means
    % Only valid if each feature is Gaussian!
    % [gmm.z, gmm.mu_mean] = kmeans(X, K);
        
end
