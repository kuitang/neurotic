function [ gmm ] = make_gmm_prior( X, K )
    [N, D] = size(X);
    % Want a STRONG pri for a HIGH variance (LOW precision).
    % So TINY gamma.
    % Got the settings from this Wikipedia page..    
    gmm = struct('K', K, ...                 
                 'prior_mean', [100 100 0.75], ...
                 'prior_cov', cov(X), ...
                 'prior_dof', 10, ...
                 'prior_scale', 10, ...
                 's_z', randsample(K, N, true), ...
                 'prior_mix', []);
    
    gmm.prior_mix = zeros(K, 1);
    gmm.prior_mix(1) = 0.4;
    gmm.prior_mix(2:end) = 0.6 / (K - 1);        
    
    %gmm.prior_mix = gmm.prior_mix * 0.0001 * N;    
    
    % Initialize with k-means
    % Only valid if each feature is Gaussian!
    % [gmm.z, gmm.mu_mean] = kmeans(X, K);    
    
    gmm.n = zeros(K, 1);    
    for n = 1:N
        k = gmm.s_z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end    
    
end
