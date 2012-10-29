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
%     gmm.s_z(n)                  ~ Discrete(mix_)
%     gmm.X(n:) | gmm.s_z(n) == k ~ N(mean_(k,:), cov_(:,:,k))   
%
% Sample log-likelihood:
%
% l(gmm.s_z | X) = sum(
    
    [N D] = size(X);
    
    % Sanity checks    
    assert(sum(gmm.n) == N);        

    % Sample statistics
    [MK, SK, NMK, NSK] = cond_moments(X, gmm.K, gmm.z);
    
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
        gmm.cov(:,:,k);
        [T p] = cholcov(gmm.cov(:,:,k));
        assert(p == 0, 'gmm.cov(:,:,k) was not PD!');
    end        
    
    %% Posterior predictive
    % Murphy (228)    
    pred_dof  = gmm.dof - D + 1;
    coef = (gmm.scale + 1) ./ (gmm.scale .* pred_dof);
    for k = 1:gmm.K        
        pred_corr(:,:,k) = coef(k) * gmm.cov(:,:,k);
    end
    
    %% Z          

    % Prepare the mvts
    
    % This MVT part is shitty... it SHOULD mix!
    % The points should NOT belong to the same class always!
    % PLOT some MVTs. Maybe your code is shit after all.
    %
    % This is why you are failing:
    % - Iteration 1: everyone has the same mean
    %   - SOMEHOW, one cluster gets all the points due to this likelihood.
    % - Iteration 2: this cluster gets all the points, so the variance is
    %   (the conjugate Bayesian update makes sense: the more points, the
    %   more you know.)
    %   - But the cluster selection is unstable. Once the variance
    %   decreases (substantially), a different cluster becomes more likely.
    %  - Iteration 3: Oscillation.
    %
    %  OBVIOUSLY THE LIKELIHOOD IS BORKED!
    %  Do you not understand MVT? The *parameters* seem to make sense. It's
    %  your likelihood that... doesn't?
    mvtparams = cell(gmm.K, 1);    
    for k = 1:gmm.K        
        mvtparams{k} = make_mvt(gmm.mean(k,:), pred_corr(:,:,k), pred_dof(k));
    end
    
    % sampling is only done here, so only these terms contribute to like
    gmm.loglike = 0;
    idxs = randperm(N);
    for nn = 1:N
        n = idxs(nn);
        k = gmm.z(n);
        % Temporarily remove ourselves
        gmm.n(k) = gmm.n(k) - 1;
        
        % Yu (4.11)
        z_prior = (gmm.n + gmm.prior_mix / gmm.K) ./ (N + gmm.prior_mix - 1);
        x_like = zeros(gmm.K, 1);
        
        % Background model
        %
        % Background is white; use a triangular distribution
        %x_like(1) = 1/N * 2 * (1 - X(n,3));                        
%         s = sqrt(diag(C));
% if (any(s~=1))
%     C = C ./ (s * s');
% end
        for k = 1:gmm.K
            % Hotspot
            x_like(k) = fast_mvtpdf(mvtparams{k}, X(n,:));            
            %x_like(k) = mvtpdf(X(n,:) - gmm.mean(k,:), pred_corr(:,:,k), pred_dof(k));                        
            %assert(abs(x_like(k) - check_like) < 1e-8, ...
            %    ['likelihoods not within 1e-8: ' num2str(x_like(k) - check_like)]);
        end                        
        
        % Yu (4.14)        
        z_pdf   = z_prior .* x_like;
        
        % Sanity check:        
        %assert(min(z_pdf) >= 10e-10*max(z_pdf));
        
        % Unnormalized inverse cdf sampling        
        z_cdf = cumsum(z_pdf);        
        k_new = find(z_cdf > z_cdf(end)*rand(1), 1);
        gmm.loglike = gmm.loglike + log(z_pdf(k_new));
        gmm.z(n) = k_new;
        gmm.n(k_new) = gmm.n(k_new) + 1;                
    end        

end
