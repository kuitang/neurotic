classdef NormalGamma < matlab.mixin.Copyable & OnlineDistribution
%NormalGamma Scalar version of NormalWishart
%
% References:
% [1] Kevin Murphy, Conjugate Bayesian analysis of the Gaussian
%     distribution, www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
    
    properties
        prior_mean, prior_n, prior_shape, prior_rate;        
        data_n
        
        post_mean, post_n, post_shape, post_rate;
        pred_mean, pred_dof, 
        
        post_mean,      post_chol, post_dof, post_n;
        
        pred_mvtparams, pred_chol, pred_dof, pred_n
    end
    
    properties (Dependent)
        pred_mean, post_cov, pred_cov
    end            
    
    methods
        function pm = get.pred_mean(o)
            pm = o.post_mean;
        end
        
        function poc = get.post_cov(o)
            poc = (o.post_chol' * o.post_chol);
        end
        
        function prc = get.pred_cov(o)
            prc = (o.pred_chol' * o.pred_chol);
        end
        
        function o = NormalWishart(prior_mean_or_rhs, prior_cov, prior_dof, prior_n)
            if nargin == 1 % copy constructor
                % REMEMBER ME!!! Standard trick to copy fields. 
                fns = properties(prior_mean_or_rhs);
                for i=1:length(fns)
                    o.(fns{i}) = rhs.(fns{i});
                end
            else                
                o.prior_mean = prior_mean_or_rhs;
                o.prior_cov  = prior_cov;
                o.prior_dof  = prior_dof;
                o.prior_n    = prior_n;

                o.dim = size(prior_cov, 1);               
                
                % IMPORTANT: Make the distribution ready for prediction and
                % update
                o.fit([]);
            end
        end                
        
        % All of the Bayesian update equations take constant time (wrt d).
        % Online updating really just updates the data covariance.
        %
        % Make sure the data covariance computation is dominating the
        % runtime before going online.
        %
        % The real kicker is are likelihood computation, methinks.
        function fit(o, X)
% fit(X) computes posterior and predictive parameters. Observations in rows
            [o.data_n, dim] = size(X);
            
            assert(o.data_n == 0 || dim == size(o.prior_cov, 1));            

            data_mean = mean(X, 1);
            data_cov  = cov(X, 1); 
            
            % Special case for empty data
            if o.data_n == 0
                data_mean = zeros(1, o.dim);
                data_cov = zeros(o.dim);
            elseif o.data_n == 1
                data_cov = zeros(o.dim);
            end                        
            
            % Unfortunately, MATLAB's method calls are too slow and we have
            % to resort to copying and pasting to eliminate much overhead
            % here
            
            % Murphy (222) to (226)
            o.post_n    = o.prior_n   + o.data_n;
            o.post_dof  = o.prior_dof + o.data_n;
            o.post_mean = 1/o.post_n * (o.prior_n * o.prior_mean + o.data_n * data_mean);
            
            % A row vector
            mean_diff = o.prior_mean - data_mean;
            post_cov  = o.prior_cov + (o.data_n*data_cov) + ...
                        (o.prior_n*o.data_n) / (o.prior_n + o.data_n) * (mean_diff'*mean_diff);
            
            % For the future rank-1 updates                        
            o.post_chol = chol(post_cov);
            
            % Murphy (228)            
            o.pred_dof = o.post_dof - o.dim + 1;            
            pred_cov = (o.post_n + 1) / (o.post_n * o.pred_dof) * post_cov;
            
            o.pred_chol = chol(pred_cov);
            o.pred_mvtparams = make_mvt(o.pred_mvtparams, o.post_mean, o.pred_chol, o.pred_dof);
        end        
        
        % Online updates. Currently only online updates for mean and cov.
        % Eventually need to update mvtparams...
        function add_point(o, x)        
            % Ross (6)
            % We mutate in place, so ORDER MATTERS!
            dx          = (o.post_mean - x)';
            c           = sqrt(o.post_n / (o.post_n + 1));                                    
            o.post_chol = cholupdate(o.post_chol, c * dx);
            
            o.post_mean = (o.post_n * o.post_mean + x) ./ (o.post_n + 1);            
            
            o.post_n    = o.post_n + 1;
            o.post_dof  = o.post_dof + 1;
            o.data_n    = o.data_n + 1;
                       
            % Murphy (228) combined with Ross (6)
            o.pred_dof  = o.pred_dof + 1;            
            
            % Note that if A = L'*L, then c^2 * A = (c*L)' * (c*L).
            % So we square-root this coefficient.
            c = (o.post_n + 1) / (o.post_n * o.pred_dof);
            o.pred_chol = sqrt(c) * o.post_chol;
            
            o.pred_mvtparams = make_mvt(o.pred_mvtparams, o.post_mean, o.pred_chol, o.pred_dof);
        end
        
        function remove_point(o, x)
            if o.data_n <= 1
                o.fit([]);
            else
                % Downdates; solve for the inverses in Ross (6)
                o.data_n    = o.data_n - 1;
                o.post_dof  = o.post_dof - 1;
                o.post_n    = o.post_n - 1;
                
                o.post_mean = ((o.post_n + 1) * o.post_mean - x) / o.post_n;                
                
                dx          = (o.post_mean - x)';
                c           = sqrt(o.post_n / (o.post_n + 1));
                o.post_chol = cholupdate(o.post_chol, c * dx, '-');                
                
%                 svec = zeros(o.dim, 1);
%                 cvec = zeros(o.dim, 1);
%                 wvec = zeros(o.dim, 1);
%                 choldnrk1({o.post_chol, [1 1 o.dim o.dim], 'U '}, ...
%                           c * dx, cvec, svec, wvec);
                
                %assertElementsAlmostEqual(naive_post_chol, o.post_chol);
                
                
                % Same Murphy's equations again                
                o.pred_dof  = o.pred_dof - 1;
                c = (o.post_n + 1) / (o.post_n * o.pred_dof);
                o.pred_chol = sqrt(c) * o.post_chol;
                
                o.pred_mvtparams = make_mvt(o.pred_mvtparams, o.post_mean, o.pred_chol, o.pred_dof);
            end
        end
        
        function [p] = pred_like(o, X)
% like(X) returns the posterior predictive likelihood of X
% TODO: Figure out the API for the non-conjugate case. Store the state as
% fields?
%
% Really need a distribution class.
            p = fast_mvtpdf(X, o.pred_mvtparams);
        end
        
        function [p] = pred_like_scalar(o, x)
            p = fast_mvtpdf_scalar(x, o.pred_mvtparams);
        end
        
        function [X] = sample_prior_pred(o, n)
% X = sample_prior(n) samples one mean and cov and n iid datapoints
            C  = iwishrnd(o.prior_cov, o.prior_dof);
            mu = mvnrnd(o.prior_mean, C);
            X  = mvnrnd(mu, C, n);
        end
        
        function [X] = sample_posterior_pred(o, n)
% X = sample_posterior(n) samples one mean and cov and n iid datapoints
            C  = iwishrnd(o.post_cov, o.post_dof);
            mu = mvnrnd(o.post_mean, C);
            X  = mvnrnd(mu, C, n);
        end
    end
    
end
