function [ mm ] = imm_iter( mm, X )
% Independent mixture model. Spatial coordinates of X (columns 1 and 2)
% full-covariance Gaussian. Remaining coordinates arbitrary.
%
% mm - mixture model structure
%      Expected fields:
%      .s_z - N x 1 vector of indicators (in 1:K)
%      .prior_mean, .prior_cov, .prior_scale, .prior_dof,
%      .mean, .cov - parameters for Gaussian part
%      prob = .background_like(mm, X) - background model. Remember to
%      normalize for spatial dimensions (1/N if you don't account for
%      space)
%      .feature_likes - cell array of function handles
%      prob = f(mm, k, X)
%
% X  - N x D matrix; first two coordinates are spatial
%
% Model
%
% mm.s_x denotes a sampled variable
% mm.p   denotes a parameter
% y_     denotes a collapsed variable
% 
% mix_ ~ Dir(mm.prior_mix)
% 
% for k=1:mm.K
%     mean_(k,:), cov_(:,:,k) ~ NW(mm.prior_mean, mm.prior_cov,
%                                  mm.prior_scale, mm.prior_dof)
%     (Parameters: mm.mean, mm.cov, mm.scale, mm.dof)
% 
% for n=1:N
%     mm.s_z(n)                   ~ Discrete(mix_)
%     mm.X(n,:) | mm.s_z(n) == k  ~ mm.background_like(mm, X)     if k == 1
%                                   N(mean_k(k,:), cov_(:,:,k)) * prod(mm.feature_likes{i}(mm, X)) 
%                                       if k >= 2
%
% Sample log-likelihood:

    [N D] = size(X);
    
    % Sanity checks    
    assert(sum(mm.n) == N);
    cov_sz = size(mm.prior_cov);
    assert(cov_sz(1) == 2);
    assert(cov_sz(2) == 2);    
    
    spatialX = X(:,1:2);
    
    % Sample statistics
    [MK SK NMK NSK] = cond_moments(spatialX, mm.K, mm.s_z);
    
    % Posterior mean
    % Bishop (10.60)
    mm.scale = mm.prior_scale + mm.n;
    % Bishop (10.61)
    mm.mean = bsxfun(@plus, mm.prior_scale .* mm.prior_mean, NMK);
    mm.mean = bsxfun(@times, 1 ./ mm.scale, mm.mean);        
    
    % Posterior covariance
    % Bishop (10.63)
    mm.dof = mm.prior_dof + mm.n;    
        
    km = (mm.prior_scale * mm.n) ./ (mm.prior_scale + mm.n);
    dm = MK - mm.mean;
    % todo: learn a tensor library
    for k = 1:mm.K
        mm.cov(:,:,k) = mm.prior_cov + NSK(:,:,k);
        mm.cov(:,:,k) = mm.cov(:,:,k) + km(k) * (dm' * dm);        
        [T p] = cholcov(mm.cov(:,:,k));
        assert(p == 0, 'mm.cov(:,:,k) was not PD!');
    end            
    
    % Posterior predictive
    % Murphy (228)    
    pred_dof  = mm.dof - D + 1;
    coef = (mm.scale + 1) ./ (mm.scale .* pred_dof);
    for k = 1:mm.K        
        pred_cov(:,:,k) = coef(k) * mm.cov(:,:,k);
    end
    
    % Data class-likelihoods       
    %
    % We will need a likelihood for each sample and each class, so
    % precompute them here.
    x_like = zeros(N, mm.K);

    % Class 1 is background (special case)
    x_like(:,1) = mm.background_like(mm, X);
    assert(all(x_like(:,1) > 0));
    
    Nfeatures = numel(mm.feature_likes);
    for k = 2:mm.K                
        % Spatial likelihood for each point and each class
        mvtparams   = make_mvt(mm.mean(k,:), pred_cov(:,:,k), pred_dof(k));
        x_like(:,k) = fast_mvtpdf(mvtparams, X(:,1:2));                       
        
        % Additional features, independent
        for nf = 1:Nfeatures
            dim = nf + 2;
            x_like(:,k) = x_like(:,k) .* mm.feature_likes{nf}(mm, k, X(:,dim));
        end        
    end        
    
    %% Sample Z
    mm.loglike = 0;
    idxs = randperm(N);
    for nn = 1:N
        n = idxs(nn);
        k = mm.s_z(n);
        % Temporarily remove ourselves
        mm.n(k) = mm.n(k) - 1;
        
        % Yu (4.11)
        z_prior = (mm.n + mm.prior_mix / mm.K) ./ (N + mm.prior_mix - 1);        

        % Yu (4.14)        
        z_pdf   = z_prior' .* x_like(n,:);
                
        % Unnormalized inverse cdf sampling        
        z_cdf = cumsum(z_pdf);        
        k_new = find(z_cdf > z_cdf(end)*rand(1), 1);                
        
        mm.loglike = mm.loglike + log(z_pdf(k_new));
        mm.s_z(n) = k_new;
        mm.n(k_new) = mm.n(k_new) + 1;                
    end        

end
