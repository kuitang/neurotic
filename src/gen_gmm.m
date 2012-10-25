function [ X, pi_, gmm ] = gen_gmm( mu, lam, N)
    % [ X, pi_, gmm ] = gen_gmm( mu, lam, N )
    [K, D] = size(mu);

    gmm = struct('z', zeros(N, 1), 'K', K, 'mu', mu, 'lam', lam);
    X = zeros(N, D);
    pi_ = rand(K, D);
    pi_ = pi_ / sum(pi_);
    pi_cdf = cumsum(pi_);
    
    for n = 1:N
        % Unnormalized inverse cdf sampling
        k = find(pi_cdf > rand(1), 1);
        gmm.z(n) = k;
        X(n, :) = normrnd(mu(k,:), 1 ./ sqrt(lam(k,:)));
    end

end

