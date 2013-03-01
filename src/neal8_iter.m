function [ mdp, state, ll ] = neal8_iter( mdp, state, params )
% Gibbs sampling iteration for Algorithm 8 in Neal 2000 [1]
%
% params is a struct with fields:
%   m - number of trial clusters (try 3)
%   b - scale for Wishart proposal
%   precision_slice_sample - slice sample?
%   mean_slice_sample   - slice sample?
%   bg_Z - normalizing constant for the background likelihood
%   graph_dist - a GraphDist loaded with 
%
% mdp.state is a n_clusters struct matrix with fields:
%   mean, prec_chol, rate - sampled posterior parameters; 4x1, 4x4, scalar.
%   bg_rate
%   
%   mdp.state(1).
%
% mdp is an MDPLight with a field prior, a struct with fields:
%   rate_shape - 
%   rate_rate  -
%   mean      -
%   n         -
%   dof       -
%   prec_chol -
%   bg_shape
%   bg_rate
%
% References:
% [1] Radford Neal, Markov Chain Sampling Methods for Dirichlet Process
%     Mixture Models
            
    ll = 0;    
    for n = 1:mdp.N
        k_old = mdp.remove_point(n);
        %assert(k_old ~= 1, 'Cannot remove the background class!');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare the auxiliary clusters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mdp.cluster_counts(k_old) == 0
            % Singleton: consider the current distribution as a new
            % candidate
            m_start = mdp.n_clusters + 2;
            m_end = m_start + params.m - 1;            
        else
            m_start = mdp.n_clusters + 1;
            m_end = m_start + params.m - 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Pseudo-sample likelihood parameters from priors  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        for k = m_start:m_end    
            [prec, mdp.state(k).prec_chol] = wishrnd([] , mdp.prior.dof, mdp.prior.prec_chol);
            %wish_loglike = wishpdfln(prec, minka_prior_dof, minka_prior_prec);
            
            mu_prec_chol = sqrt(mdp.prior.n) * mdp.state(k).prec_chol;
            mdp.state(k).mean = randnorm(1, mdp.prior.mean, mu_prec_chol);
            %norm_loglike = mvnormpdfln(mdp.state(k).mean, mdp.prior.mean, mu_prec_chol, 'inv');
            
            mdp.state(k).rate = gamrnd(mdp.prior.rate_shape, 1 / mdp.prior.rate_rate);
            %exp_loglike = -gamlike([mdp.prior.rate_shape 1 / mdp.prior.rate_rate], mdp.state(k).rate);
        end
         
        % Actually, we don't need the log likelihood. These were just
        % pseudo-samples. We only add in the log likelihood of the discrete
        % z we sample.
        %
        %ll = ll + wish_loglike + norm_loglike + exp_loglike;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample z_n
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Nclusters = length(mdp.cluster_counts) + params.m;
        z_pdf = [mdp.cluster_counts ; mdp.concentration ./ (params.m * ones(params.m,1) )] ...
                ./ ( mdp.N - 1 + mdp.concentration );
        
        x_pdf = zeros(Nclusters, 1);
        
        % Calculate the background
        bg_exp_like = gampdf(mdp.X(n,3), mdp.prior.bg_shape, 1 / mdp.state(1).bg_rate);
        x_pdf(1) = bg_exp_like / params.bg_Z;
        
        for k = 2:params.m            
            s = mdp.state(k);
            
            % NO ;;
            norm_like = mvnormpdf(mdp.X(n,1:4)', s.mean, s.prec_chol, 'inv');            
            exp_like  = exppdf(params.graph_dist.dist(mdp.X(n,:), s.mean), 1 / s.rate);
            x_pdf(k) = norm_like * exp_like;
        end
            
        % NO ;;
        pdf = z_pdf .* x_pdf;
        
        % Unnormalized inverse cdf sampling        
        cdf = cumsum(pdf);        
        k_new = find(cdf > cdf(end)*rand(1), 1);
        mdp.assign(n, k_new);
        mdp.prune_clusters();
                
        ll = ll + log(pdf(k_new));       
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update cluster parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Translate Murphy/MATLAB/Bishop parameterization to Minka's
    % see Murphy (291).
    
    % We use Minka's wishpdfln to use for the sampling distribution.
    % While the Wishart distribution takes a template precision and
    % spits out another precision (the higher the degrees of freedom,
    % the more precise the precision).
    %
    % Thus, notationally, most (including MATLAB's wishpdf)
    % parameterizations of the Wishart distribution take a template
    % precision prec. Minka's precision takes a parameter 
    %
    % B = -0.5*inv(prec).    
    % a = nu / 2
    %
    % Minka expects an inverted parameter. But it's easy to invert Cholesky
    % factors.
    
    inv_prec_chol = inv(mdp.prior.prec_chol);
    minka_prior_prec = 0.5 * (inv_prec_chol * inv_prec_chol');
    minka_prior_dof  = mdp.prior.dof / 2;    
                  
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample background posteriors (conjugate)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bgidxs = mdp.cluster_assigns == 1;
    post_bg_shape = mdp.prior.bg_shape + sum(bgidxs);
    post_bg_rate  = mdp.prior.bg_rate + mdp.X(bgidxs,3);
    
    mdp.state(1).bg_rate = gamrnd(post_bg_shape, 1 / post_bg_rate);
    ll = ll - gamlike([post_bg_shape,  1 / post_bg_rate], mdp.state(1).bg_rate);
    
    for k = 2:mdp.n_clusters        
        kidxs = mdp.cluster_assigns == k;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample rate (conjugate)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dists = params.graph_dist.dist(mdp.X(kidxs,:), s.mean);        
        post_rate_shape = mdp.prior.rate_shape + length(dists);
        post_rate_rate  = mdp.prior.rate_rate  + sum(dists);
        
        mdp.state(k).bg_rate = gamrnd(post_rate_shape, 1 / post_rate_rate);
        ll = ll - gamlike([post_rate_shape, 1 / post_rate_rate], mdp.state(k).rate);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample precision (MH)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        prec_chol_loglike = @(prec_chol) log_xlike(mdp.X(kidxs,1:4)', s.mean, prec_chol, mdp) + ...
                            wishpdfln(prec_chol' * prec_chol, minka_prior_prec, minka_prior_dof);
        prec_chol_logprop = @(new_chol, prev_chol) wishpdfln(new_chol' * new_chol, params.b * prev_chol, minka_prior_dof);
        prec_chol_rnd = make_prec_chol_rnd();

        if params.precision_slice_sample
            % TODO: READ THE DAMN PAPER AND FIGURE OUT THE PARAM SETTINGS!
            error('Not implemented');
        else
            % To keep the Markov chain invariant, we start MH at our
            % current sample. But our functions take Cholesky factors.            
            [prec_chol_samp, accept_rate] = ...
                mhsample(s.prec_chol, 1, 'logpdf', prec_chol_loglike, ...
                         'logproppdf', prec_chol_logprop, 'proprnd', prec_chol_rnd);
            
            % NO ;;
            accept_rate
              
            mdp.state(k).prec = prec_chol_samp' * prec_chol_samp;            
        end        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample mean (MH)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mu_loglike = @(mu) log_xlike(mdp.X(kidxs,:), mu, s.prec_chol, mdp) + ...
                           mvnormpdfln(mu, mdp.prior.mean, prior_mean_prec_chol, 'inv');
        mu_logprop = @(new_mu, prev_mu) mvnormpdfln(new_mu, prev_mu, params.mu_prop_chol, 'inv');
        mu_rnd     = @(prev_mu) randnorm(1, prev_mu, params.mu_prop_chol);                        
        
                
        if params.mean_slice_sample
            error('Not implemented');
        else
            [mu_samp, accept_rate] = ...
                mhsample(s.mean', 1, 'logpdf', mu_loglike, ...
                         'logproppdf', mu_logprop', 'proprnd', mu_rnd);
            
            % NO ;;
            accept_rate            
            mdp.state(k).mean = mu_samp;
        end
    end

end

function rn = make_prec_chol_rnd()
    function c = inner(prev_chol)
        [~, c] = wishrnd([], mdp.prior_mean_prec.prior_dof, params.b * prev_chol);
    end

    rn = @inner;
end

function ll = log_xlike(X, mu, prec_chol, mdp)
% ll = log_xlike(X, mu, prec)
%
% X    - DxNk data matrix
% mu   - Dx1     COLUMN vector (mean)
% prec - DxD     FULL matrix, not just Cholesky factor (precision)
%
% prec is the output of the Wishart prior. It is NOT scaled by kappa = 
% prior_n. This routine does the scaling inside.
%
% Log likelihood for class k. Uses Minka's mvnormpdf with FULL precision or
% covariance (precision if do_inv is 'inv'). NOT using cholesky factors.

    N = size(X, 2);    
    scaled_prec_chol = sqrt(mdp.prior.n) * prec_chol;        
    ll = sum(mvnormpdfln(X, mu, scaled_prec_chol, 'inv'));
end
