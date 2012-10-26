function [ gmm ] = gmm_gibbs_iter( gmm, X )
    % Full covariance collapsed Gibbs sampling    
    [N D] = size(X);
    
    % Sanity checks    
    assert(sum(gmm.n) == N);        
    
    % This Normal-Wishart thing is annoying; let's just do the normal and
    % the Wishart part separately.
    
    %% LAMBDA
    [MK, SK] = cond_moments(X, gmm.z);
    % Murphy (224)    
    for k = 1:gmm.K
        SK(:,:,k) = gmm.n(k) * SK(:,:,k);
    end
    % Murphy (225)
    lam_dof = gmm.lam_dof + gmm.n; 
        
    shape_scale = (gmm.scale_prior * gmm.n) ./ (gmm.scale_prior + gmm.n);
    for k = 1:gmm.K
        dm = (gmm.mu_mean(k,:) - MK(k,:));
        assert(all(size(dm) == [1 D]));
        
        lam_shape = gmm.lam_shape + SK(:,:,k) + shape_scale(k) * (dm'*dm)
        
        gmm.lam(:,:,k) = wishrnd(lam_shape, lam_dof(k));        
    end                    
    
    %% MU
    % Murphy (226)
    gmm.scale = gmm.scale_prior + gmm.n;
    % Murphy (222)        
    mu_mean_top = (gmm.scale_prior * gmm.mu_mean + bsxfun(@times, gmm.n, MK));
    mu_mean_bot = (gmm.scale_prior + gmm.n);
    mu_mean     = bsxfun(@times, mu_mean_top, 1 ./ mu_mean_bot);        
    
    for k = 1:gmm.K        
        gmm.mu(k,:) = mvnrnd(mu_mean(k,:), inv(gmm.scale(k) * gmm.lam(:,:,k)));
    end
    
    gmm.mu
   
    %% Z
    
    % pi was integrated out    
    % TODO: SPECIAL CASE THE BACKGROUND MODEL
    
    % Upper triangle cholesky, for faster mvnormpdf
    inv_stdev = zeros(D, D, gmm.K);
    log_Z = zeros(gmm.K, 1);
    for k = 1:gmm.K
        inv_stdev(:,:,k) = chol(gmm.lam(:,:,k));
        log_Z(k) = sum(log(diag(inv_stdev(:,:,k)))) - 0.5 * D * log(2 * pi);
    end
    
    gmm.mu
    
    for n = 1:N
        k = gmm.z(n);
        % Temporarily remove ourselves
        gmm.n(k) = gmm.n(k) - 1;
        
        % Yu (4.11)
        z_prior = (gmm.n + gmm.alpha / gmm.K) ./ (N + gmm.alpha - 1);
        x_like = zeros(gmm.K, 1);
        for k = 1:gmm.K
            % Hotspot -- copied from lightspeed            
            check_like = mvnpdf(X(n,:), gmm.mu(k,:), inv(gmm.lam(:,:,k)));
            
            dx = X(n,:) - gmm.mu(k,:);
            dx_std = dx * inv_stdev(:,:,k);
            log_like = log_Z(k) - 0.5 * (dx_std*dx_std');
            x_like(k) = exp(log_like);
            assert(x_like(k) > 0);
            [x_like(k), check_like];
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
