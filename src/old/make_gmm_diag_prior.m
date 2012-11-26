function [ gmm ] = make_gmm_diag_prior( X, K )
    [N, D] = size(X);
    % Want a STRONG prior for a HIGH variance (LOW precision).
    % So TINY gamma.
    % Got the settings from this Wikipedia page..    
    gmm = struct('K', K, ...
                 'mu_prec', repmat([1e-4 1e-4 1], K, 1), ... % Means are robust from the data.
                 'mu_mean', repmat([100 100 0.75], K, 1), ...
                 'lam_shape', repmat([4 4 6400], K, 1), ... % See board
                 'lam_scale', repmat([1000 1000 100], K, 1), ...
                 'z', randsample(K, N, true), ...       % This doesn't matter
                 'alpha', []);        % Each component starts 20% full. Essential to prevent collapsing to 0. (That's not true)
    % Gibbs sampling needs to start with something...
    gmm.lam = gamrnd(gmm.lam_shape, 1 ./ gmm.lam_scale);
    assert(all(size(gmm.lam) == [K D]));        
    
    gmm.alpha = zeros(K, 1);
    gmm.alpha(1) = 0.4
    gmm.alpha(2:end) = 0.6 / (K - 1);
    gmm.alpha = gmm.alpha * 0.001 * N;
    
    % Initialize with k-means
    % Only valid if each feature is Gaussian!
    %[gmm.z, gmm.mu_mean] = kmeans(X, K);    
    
    gmm.n = zeros(K, 1);    
    for n = 1:N
        k = gmm.z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end    
    
end

