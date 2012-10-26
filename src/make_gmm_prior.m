function [ gmm ] = make_gmm_prior( X, K )
    [N, D] = size(X);
    % Want a STRONG prior for a HIGH variance (LOW precision).
    % So TINY gamma.
    % Got the settings from this Wikipedia page..    
    gmm = struct('K', K, ...
                 'lam_shape', [], ...
                 'lam', [], ...
                 'lam_dof', 20, ...                      % Upper triangle
                 'mu_mean', rand(K, D), ...             % Means are robust from the data.
                 'scale_prior', 0.1, ...                % We never want much precision. Our Bayesian update mayes Wisharts very precise.
                 'z', randsample(K, N, true), ...       % This doesn't matter
                 'alpha', 0.2*(N/K)*ones(K, 1));        % Each component starts 20% full. Essential to prevent collapsing to 0.
                 
    % This is all I know
    gmm.lam_shape =  1e-15 * [ 1 1 0.3
                      1 1 0.3
                      0.3 0.3 0.1 ];
                  
    % Gibbs sampling needs to start with something...
    gmm.lam = zeros(D, D, K);
    for k = 1:K        
        gmm.lam(:,:,k) = wishrnd(gmm.lam_shape, gmm.lam_dof);
    end

    
    % Initialize with k-means
    % Only valid if each feature is Gaussian!
    [gmm.z, gmm.mu_mean] = kmeans(X, K);    
    
    gmm.n = zeros(K, 1);    
    for n = 1:N
        k = gmm.z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end    
    
end

