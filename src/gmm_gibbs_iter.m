function [ gmm ] = gmm_gibbs_iter( gmm, X )
    % Full covariance collapsed Gibbs sampling    
    [N D] = size(X);
    
    % Sanity checks    
    assert(sum(gmm.n) == N);        
            
    %% ?
    [MK, SK] = cond_moments(X, gmm.z);    
    % Murphy (225)
    lam_dof = gmm.lam_dof + gmm.n;         
            
    shape_scale = (gmm.kappa * gmm.n) ./ (gmm.kappa + gmm.n);
    for k = 1:gmm.K
        dm = (gmm.mu_mean(k,:) - MK(k,:));
        assert(all(size(dm) == [1 D]));
        
        gmm.lam_cov
        disp('C')
        SK(:,:,k) * gmm.n(k)
        disp('dm')
        shape_scale(k) * (dm'*dm)
        
        lam_shape = gmm.lam_cov + SK(:,:,k) + shape_scale(k) * (dm'*dm);
        
        % Murphy (223)
        gmm.lam(:,:,k) = wishrnd(inv(lam_shape), lam_dof(k)) / lam_dof(k);
    end    
    
    %% ?
    % Murphy (226)
    gmm.scale = gmm.kappa + gmm.n;
    % Murphy (222)        
    mu_mean_top = (gmm.kappa * gmm.mu_mean + bsxfun(@times, gmm.n, MK));
    mu_mean_bot = (gmm.kappa + gmm.n);
    mu_mean     = bsxfun(@times, mu_mean_top, 1 ./ mu_mean_bot);        
    
    MU_VAR = inv(gmm.scale(k) * gmm.lam(:,:,k))
    
    for k = 1:gmm.K        
        gmm.mu(k,:) = mvnrnd(mu_mean(k,:), inv(gmm.scale(k) * gmm.lam(:,:,k)));
    end        
   
    %% Z
    
    % pi was integrated out    
    % TODO: SPECIAL CASE THE BACKGROUND MODEL
    
    % Upper triangle cholesky, for faster mvnormpdf
    inv_stdev = zeros(D, D, gmm.K);
    log_Z = zeros(gmm.K, 1);
    for k = 1:gmm.K
        inv_stdev(:,:,k) = chol(gmm.lam(:,:,k));
        log_Z(k) = sum(log(diag(inv_stdev(:,:,k)))) - 0.5 * D * log(2 * pi);
        CLASS_VAR = inv(gmm.lam(:,:,k))
    end        
    
    for n = 1:N
        k = gmm.z(n);
        % Temporarily remove ourselves
        gmm.n(k) = gmm.n(k) - 1;
        
        % Yu (4.11)
        z_prior = (gmm.n + gmm.alpha / gmm.K) ./ (N + gmm.alpha - 1);
        x_like = zeros(gmm.K, 1);
        for k = 1:gmm.K
            % Hotspot -- copied from lightspeed            
            %check_like = mvnpdf(X(n,:), gmm.mu(k,:), inv(gmm.lam(:,:,k)));
            
            dx = X(n,:) - gmm.mu(k,:);
            dx_std = dx * inv_stdev(:,:,k);
            log_like = log_Z(k) - 0.5 * (dx_std*dx_std');
            x_like(k) = exp(log_like);
            %assert(x_like(k) > 0);
            %assert(abs(x_like(k) - check_like) < 1e-8, x_like(k) - check_like);
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
