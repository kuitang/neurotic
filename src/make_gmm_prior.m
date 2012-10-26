function [ gmm ] = make_gmm_prior( X, K )
    [N, D] = size(X);
    % Want a STRONG prior for a HIGH variance (LOW precision).
    % So TINY gamma.
    % Got the settings from this Wikipedia page..    
    gmm = struct('K', K, ...
                 'lam_cov', [], ...
                 'lam', [], ...
                 'lam_dof', 10, ...                      % Upper triangle
                 'mu_mean', rand(K, D), ...             % Means are robust from the data.
                 'kappa', 0.1*(N/K), ...                
                 'z', randsample(K, N, true), ...       % This doesn't matter
                 'alpha', 0.2*(N/K)*ones(K, 1));        % Each component starts 20% full. Essential to prevent collapsing to 0.
                 
    % This is all I know
    gmm.lam_cov =  100 * [ 1 -0.1 0.2; -0.1 1 0.2; 0.2 0.2 0.3];
                  
    % Gibbs sampling needs to start with something...
    gmm.lam = zeros(D, D, K);
    for k = 1:K        
        gmm.lam(:,:,k) = wishrnd(inv(gmm.lam_cov), gmm.lam_dof) / gmm.lam_dof;
    end
    
    gmm.lam

    
    % Initialize with k-means
    % Only valid if each feature is Gaussian!
    [gmm.z, gmm.mu_mean] = kmeans(X, K);    
    
    gmm.n = zeros(K, 1);    
    for n = 1:N
        k = gmm.z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end    
    
end

