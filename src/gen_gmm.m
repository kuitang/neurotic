function [ X, pi, gmm ] = gen_gmm( mu, lam, N )
    
    [K, D] = size(mu);

    gmm = struct('z', zeros(N, 1), 'K', K, 'mu', mu, 'lam', lam);
    X = zeros(N, 2);
    pi = rand(K, 1);    
    pi_cdf = cumsum(pi)    
    
    for n = 1:N
        % Unnormalized inverse cdf sampling
        k = find(pi_cdf > pi_cdf(end)*rand(1), 1)
        gmm.z(n) = k;
        X(n, :) = normrnd(mu(k,:), 1 ./ sqrt(lam(k,:)));
    end

end

