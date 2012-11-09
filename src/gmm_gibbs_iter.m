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
    %assert(sum(gmm.n) == N);        

    % Class 1 is background.
    gmm.x_like(:,1) = 1/gmm.N * gmm.background_pdf(X(:,3));
    %assert(all(x_like(:,1) > 0));    
    
    %% Sample
    gmm.loglike = 0;
    idxs = randperm(N);
    for nn = 1:N
        n = idxs(nn);
        k = gmm.s_z(n);                

        if k > 0        
            % If we belong to a cluster, temporarily remove ourselves
            gmm.s_z(n) = 0;
            gmm.n(k)   = gmm.n(k) - 1;
            if gmm.n(k) > 0
                gmm = gmm_recompute_cluster(gmm, X, k);
            else
                % Radford Neal says to keep the old parameters. To
                % implement later.                
                gmm.empty_clusters = [gmm.empty_clusters k];    
                gmm.K = gmm.K - 1;
            end            
        end
        
        nzidxs = gmm.n > 0;
        assert(sum(nzidxs) == gmm.K, 'gmm.K does not track # of nonzero classes');
        
        % Neal (3.7); Algorithm 3        
        
        % Compute likelihood for the Dirichlet part. The gmm.n accounts for
        % the existing classes and the gmm.prior_conc scalar accounts for
        % the inchoate classes.
        z_prior = [ gmm.n(nzidxs) ; gmm.prior_conc ] ./ (sum(gmm.n) - 1 + gmm.prior_conc);
        
        % Compute the posterior predictive likelihoods with ourselves
        % removed, for each of the existing classes.        
        %[x_like, gmm] = gmm_mvnw_posterior_pred(gmm, X, n);
        
        % We've precomputed our new-class likelihoods, so augment here.                
        % Combine the Dirichlet and Gaussian parts.
        x_like = [ gmm.x_like(n,nzidxs) gmm.new_like(n) ];
        z_pdf  = z_prior' .* x_like;        
                
        % Unnormalized inverse cdf sampling        
        z_cdf = cumsum(z_pdf);        
        ik_new = find(z_cdf > z_cdf(end)*rand(1), 1);                      
        gmm.loglike = gmm.loglike + log(z_pdf(ik_new));
        
        % Recover the true class index
        if ik_new <= gmm.K % if we picked an existing cluster
            knz = find(nzidxs);
            k_new = knz(ik_new);            
        else
            assert(~isempty(gmm.empty_clusters), 'no free clusters left!');
            % Assign k_new to the next empty cluster
            k_new = gmm.empty_clusters(1);
            gmm.empty_clusters(1) = [];            
            gmm.K = gmm.K + 1;
        end
        
        gmm.s_z(n) = k_new;
        gmm.n(k_new) = gmm.n(k_new) + 1;         
                        
        % Officially add ourselves to the cluster
        %gmm = add_to_cluster(gmm, X, n, k_new);
        gmm = gmm_recompute_cluster(gmm, X, k_new);
        
    end
end
