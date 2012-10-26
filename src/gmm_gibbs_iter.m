function [ gmm ] = gmm_gibbs_iter( gmm, X )
    % Diagonal covariance collapsed Gibbs sampling    
    [N D] = size(X);
    
    % Sanity checks    
    assert(sum(gmm.n) == N);        
    
    %% MU
    % Murphy (29)
    n_lam = bsxfun(@times, gmm.n, gmm.lam);
    mu_prec = gmm.mu_prec + n_lam;    
    w = bsxfun(@times, gmm.n, gmm.lam ./ mu_prec);
    assert(all(size(w) == [gmm.K, D]));
    
    MX = cond_mean(X, gmm.z);
    % Murphy (30)
    mu_mean = bsxfun(@times, w, MX) + bsxfun(@times, 1 - w, gmm.mu_mean);
    assert(all(size(mu_mean) == [gmm.K, D]));
    
    gmm.mu = normrnd(mu_mean, 1 ./ sqrt(mu_prec));    
    assert(all(size(gmm.mu) == [gmm.K, D]));
    assert(all(all(~isnan(gmm.mu))));
    
    %% LAMBDA
    % Murphy (117)
    lam_shape = bsxfun(@plus, gmm.lam_shape, gmm.n/2); % Remember to adjust per-class counts
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
    D2_log_2pi = D/2 * log(2*pi);
    sum_log_lam = 1/2 * sum(log(gmm.lam), 2);
    assert(all(size(sum_log_lam) == [gmm.K 1]));
    
    for n = 1:N
        k = gmm.z(n);
        % Temporarily remove ourselves
        gmm.n(k) = gmm.n(k) - 1;
        
        % Yu (4.11)
        z_prior = (gmm.n + gmm.alpha / gmm.K) ./ (N + gmm.alpha - 1);
        x_like = zeros(gmm.K, 1);
        for k = 1:gmm.K
            % Hotspot
            X_diff = X(n,:) - gmm.mu(k,:);
            quadform = -1/2 * sum(gmm.lam(k,:) .* (X_diff.^2));
            %x_like(k) = sqrt_2pi_to_D * sqrt(prod(gmm.lam(k,:))) * exp(quadform);                        
            log_like = -D2_log_2pi + sum_log_lam(k) + quadform;
            x_like(k) = exp(log_like);

            %old_norm_like = prod(normpdf(X(n,:), gmm.mu(k,:), 1 ./ sqrt(gmm.lam(k,:))));                
            %[x_like(k), old_norm_like]
            %assert(abs(x_like(k) - old_norm_like) < 10e-5);
            
            %[k x_like(k); X(n,:); gmm.mu(k,:); 1 ./ gmm.lam(k,:)]
        end
        %x_like
        
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
