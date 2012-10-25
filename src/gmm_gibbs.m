function [ gmm ] = gmm_gibbs( X, K, niters )

    [N, D] = size(X);
    gmm = struct('K', K, 'mu_prec', ones(K, D), ...
                 'lam_shape', ones(K, D), 'lam_scale', ones(K, D), ...
                 'alpha', ones(K, 1));
    % Gibbs sampling needs to start with something...
    gmm.lam = gamrnd(gmm.lam_shape, 1 ./ gmm.lam_scale);    
    gmm.n = zeros(K, 1);    
    
    % Initialize with k-means
    %[gmm.z, gmm.mu_mean] = kmeans(X, K);    
    gmm.z = randsample(K, N, true);
    gmm.mu_mean = rand(2,2);
        
    for n = 1:N
        k = gmm.z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end
    
    for i = 1:niters
        gmm = gmm_gibbs_iter(gmm, X);
    end

end

