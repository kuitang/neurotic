function [ gmm ] = gmm_gibbs( X, K, niters )

    [N, D] = size(X);
    % Weak prior
    gmm = struct('K', K, ...
                 'mu_prec', 0.1*ones(K, D), ...
                 'mu_mean', rand(K, D), ...
                 'lam_shape', 0.1*ones(K, D), ...
                 'lam_scale', 0.1*ones(K, D), ...
                 'z', randsample(K, N, true), ...
                 'alpha', 10*ones(K, 1));
    % Gibbs sampling needs to start with something...
    gmm.lam = gamrnd(gmm.lam_shape, 1 ./ gmm.lam_scale);
    assert(all(size(gmm.lam) == [K D]));    
    
    % Initialize with k-means
    [gmm.z, gmm.mu_mean] = kmeans(X, K);    
    
    gmm.n = zeros(K, 1);    
    for n = 1:N
        k = gmm.z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end
    
    for i = 1:niters
        tic
        z_old = gmm.z;
        gmm = gmm_gibbs_iter(gmm, X);
        change = mean(abs(z_old ~= gmm.z));
        disp(['Iteration ' num2str(i) ' took ' num2str(toc) ' seconds; ' num2str(change) ' changed classes']);
        gmm.n        
    end

end

