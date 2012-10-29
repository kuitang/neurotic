function [ X, pi_, gmm ] = gen_gmm(mu, N)
    % Generate a mixture of Gaussians to test algos; mu K x D, lam D x D x K
    [K, D] = size(mu);
    
    C = zeros(D, D, K);    
    for k = 1:K
        % The stupidest way to generate covariances: sample covariances.
        % scale by D
        while true
            C(:,:,k) = D^2 * cov(rand(D,D));
            [sd t] = cholcov(C(:,:,k));
            if t == 0 && rank(sd) == D
                break;
            end
        end
    end
    
    gmm = struct('z', zeros(N, 1), 'K', K, 'n', zeros(K, 1));
    X = zeros(N, D);
    pi_ = rand(K, 1);
    pi_ = pi_ / sum(pi_);
    
    gmm.z = sample(pi_, N)';
        
    for k = 1:K
        idxs = gmm.z == k;
        gmm.n(k) = sum(idxs);
        
        [stdev, t] = cholcov(C(:,:,k));
        assert(all(size(stdev) == [D D]));
        assert(t == 0, 'not psd!');
        X(idxs,:) = randnorm(gmm.n(k), mu(k,:)', stdev)';
    end
    
    % Now construct an uninformative prior
    gmm.prior_mean  = rand(K, D);
    gmm.prior_scale = 0.1;
    gmm.prior_dof   = D + 1;
    gmm.prior_cov   = zeros(D, D, K);
    for k = 1:K
        gmm.prior_cov(:,:,k) = eye(D);
    end
    gmm.prior_mix   = 2 * ones(K, 1);

end
