function [ gmm ] = make_gmm_prior( X, K )
    [N, D] = size(X);
    % Want a STRONG prior for a HIGH variance (LOW precision).
    % So TINY gamma.
    % Got the settings from this Wikipedia page..    
    gmm = struct('K', K, ...
                 'mu_prec', 0.1*ones(K, D), ... % Means are robust from the data.
                 'mu_mean', rand(K, D), ...
                 'lam_shape', K/N*0.001*ones(K, D), ... % The density will hug 0.0001.
                 'lam_scale', 100*ones(K, D), ...       % And our variance estimates shall always be high
                 'z', randsample(K, N, true), ...       % This doesn't matter
                 'alpha', 0.2*(N/K)*ones(K, 1));        % Each component starts 20% full. Essential to prevent collapsing to 0.
    % Gibbs sampling needs to start with something...
    gmm.lam = gamrnd(gmm.lam_shape, 1 ./ gmm.lam_scale);
    assert(all(size(gmm.lam) == [K D]));
    
    % Initialize with k-means
    % Only valid if each feature is Gaussian!
    %[gmm.z, gmm.mu_mean] = kmeans(X, K);    
    
    gmm.n = zeros(K, 1);    
    for n = 1:N
        k = gmm.z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end    
    
end

