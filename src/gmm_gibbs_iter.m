function [ gmm ] = gmm_gibbs_iter( gmm, X )
    [N D] = size(X);            
    
    % Diagonal covariance
    
    %% MU        
    % Murphy (29)    
    mu_prec = gmm.mu_mean + N*gmm.lam;
    w = N*gmm.lam ./ mu_prec;
    assert(all(size(w) == [gmm.K, D]));
    
    MX = cond_mean(X, gmm.z);
    % Murphy (30)
    mu_mean = bsxfun(@times, w, MX) + bsxfun(@times, 1 - w, gmm.mu_mean);
    assert(all(size(mu_mean) == [gmm.K, D]));
    
    gmm.mu = normrnd(mu_mean, 1 ./ sqrt(mu_prec));
    gmm.mu;
    assert(all(size(gmm.mu) == [gmm.K, D]));
    assert(all(all(~isnan(gmm.mu))));
    
    %% LAMBDA
    % Murphy (117)
    lam_shape = gmm.lam_shape + N/2;        
    lam_scale = gmm.lam_scale;
        
    for n = 1:N
        k = gmm.z(n);
        X_diff = X(n, :) - gmm.mu(k, :);
        dev = 0.5 * X_diff.^2;
        %[X(n, :); gmm.mu(k, :)]
        
        % Murphy (118)
        lam_scale(k, :) = lam_scale(k, :) + dev;
    end        
    
    gmm.lam = gamrnd(lam_shape, 1 ./ lam_scale);
    assert(all(size(gmm.lam) == [gmm.K, D]));
    
    %% Z
    
    % pi was integrated out    
    % TODO: SPECIAL CASE THE BACKGROUND MODEL
    for n = 1:N
        k = gmm.z(n);
        % Temporarily remove ourselves
        gmm.n(k) = gmm.n(k) - 1;
        
        % Yu (4.11)
        z_prior = (gmm.n + gmm.alpha / k) ./ (N + gmm.alpha - 1);
        x_like = ones(gmm.K, 1);
        for k = 1:gmm.K            
            x_like(k) = mvnpdf(X(n,:), gmm.mu(k,:), 1 ./ gmm.lam(k,:));
            %[k x_like(k); X(n,:); gmm.mu(k,:); 1 ./ gmm.lam(k,:)]
        end
        
        % Yu (4.14)
        z_pdf   = z_prior .* x_like;
        % Sanity check:        
        %assert(min(z_pdf) >= 10e-10*max(z_pdf));
        
        % Unnormalized inverse cdf sampling
        z_cdf = cumsum(z_pdf);
        k_new = find(z_cdf > z_cdf(end)*rand(1), 1);        
        gmm.z(n) = k_new;
        gmm.n(k_new) = gmm.n(k_new) + 1;
    end        

end
