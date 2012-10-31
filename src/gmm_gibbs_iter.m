function [ gmm ] = gmm_gibbs_iter( gmm, X )
% Full covariance collapsed Gibbs sampling
%
% gmm.s_x denotes a sampled variable
% gmm.p   denotes a parameter
% y_      denotes a collapsed variable
% 
% mix_ ~ Dir(gmm.prior_mix)
% 
% for k=1:gmm.K
%     mean_(k,:), cov_(:,:,k) ~ NW(gmm.prior_mean, gmm.prior_cov,
%                                  gmm.prior_scale, gmm.prior_dof)
%     (Parameters: gmm.mean, gmm.cov, gmm.scale, gmm.dof)
% 
% for n=1:N
%     gmm.s_z(n)                   ~ Discrete(mix_)
%     gmm.X(n,:) | gmm.s_z(n) == k ~ 1/N * 2 * 2*X(n,:)           if k == 1
%                                    N(mean_(k,:), cov_(:,:,k))   if k >= 1
%
% Sample log-likelihood:
%

    [N D] = size(X);
    
    % Sanity checks    
    assert(sum(gmm.n) == N);        

    % Sample statistics
    [MK, SK, NMK, NSK] = cond_moments(X, gmm.K, gmm.s_z);
    
    %% MU            
    % Bishop (10.60)
    gmm.scale = gmm.prior_scale + gmm.n;
    % Bishop (10.61)
    gmm.mean = bsxfun(@plus, gmm.prior_scale .* gmm.prior_mean, NMK);
    gmm.mean = bsxfun(@times, 1 ./ gmm.scale, gmm.mean);        
    
    %% LAMBDA
    % Bishop (10.63)
    gmm.dof = gmm.prior_dof + gmm.n;    
        
    km = (gmm.prior_scale * gmm.n) ./ (gmm.prior_scale + gmm.n);
    dm = MK - gmm.mean;
    % todo: learn a tensor library
    for k = 1:gmm.K
        gmm.cov(:,:,k) = gmm.prior_cov + NSK(:,:,k);
        gmm.cov(:,:,k) = gmm.cov(:,:,k) + km(k) * (dm' * dm);        
        [T p] = cholcov(gmm.cov(:,:,k));
        assert(p == 0, 'gmm.cov(:,:,k) was not PD!');
    end            
    
    %% Posterior predictive
    % Murphy (228)    
    pred_dof  = gmm.dof - D + 1;
    coef = (gmm.scale + 1) ./ (gmm.scale .* pred_dof);
    for k = 1:gmm.K        
        pred_cov(:,:,k) = coef(k) * gmm.cov(:,:,k);
    end
    
    %% Z          
    
    % We will need a likelihood for each sample and each class, so
    % precompute them here.
    x_like = zeros(N, gmm.K);

    % Class 1 is background.     
    x_like(:,1) = gmm.background_like(gmm, X);
    assert(all(x_like(:,1) > 0));
    % Precompute likelihood for every point for every class!!!
    % That means you don't check idxs = gmm.s_z == k; idiot.
    
    for k = 2:gmm.K                
        mvtparams   = make_mvt(gmm.mean(k,:), pred_cov(:,:,k), pred_dof(k));
        x_like(:,k) = fast_mvtpdf(mvtparams, X);                       
    end        
    
    % Now, sample
    gmm.loglike = 0;
    idxs = randperm(N);
    for nn = 1:N
        n = idxs(nn);
        k = gmm.s_z(n);
        % Temporarily remove ourselves
        gmm.n(k) = gmm.n(k) - 1;
        
        % Yu (4.11)
        z_prior = (gmm.n + gmm.prior_mix / gmm.K) ./ (N + gmm.prior_mix - 1);        

        % Yu (4.14)        
        z_pdf   = z_prior' .* x_like(n,:);
                
        % Unnormalized inverse cdf sampling        
        z_cdf = cumsum(z_pdf);        
        k_new = find(z_cdf > z_cdf(end)*rand(1), 1);                
        
        gmm.loglike = gmm.loglike + log(z_pdf(k_new));
        gmm.s_z(n) = k_new;
        gmm.n(k_new) = gmm.n(k_new) + 1;                
    end        

end
