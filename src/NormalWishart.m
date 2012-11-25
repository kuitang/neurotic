classdef NormalWishart < matlab.mixin.Copyable
%NormalWishart Deterministically fit a Gaussian likelihood to its conjugate prior
%
% References:
% [1] Kevin Murphy, Conjugate Bayesian analysis of the Gaussian
%     distribution, www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
% [2] Andrew Gelman et al. Bayesian Data Analysis, Second Edition.
    
    properties
        prior_mean, prior_cov, prior_dof, prior_n, dim
        data_mean, data_cov, data_n
        post_mean, post_cov, post_dof, post_n
        pred_mean, pred_cov, pred_dof, pred_n, pred_mvtparams        
    end
    
    methods
        function o = NormalWishart(prior_mean, prior_cov, prior_dof, prior_n)
            o.prior_mean = prior_mean;
            o.prior_cov  = prior_cov;
            o.prior_dof  = prior_dof;
            o.prior_n    = prior_n;
            
            o.dim = size(prior_cov, 1);
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

            o.data_mean = mean(X, 1);
            o.data_cov  = cov(X, 1); 
            
            % Special case for empty data
            if o.data_n == 0
                o.data_mean = 0;
                o.data_cov = 0;
            end
            
            % Murphy (222) to (226)
            o.post_n    = o.prior_n   + o.data_n;
            o.post_dof  = o.prior_dof + o.data_n;
            o.post_mean = 1/o.post_n * (o.prior_n * o.prior_mean + o.data_n * o.data_mean);
            
            % A row vector
            mean_diff = o.prior_mean - o.data_mean;
            o.post_cov  = o.prior_cov + (o.data_n*o.data_cov) + ...
                        (o.prior_n*o.data_n) / (o.prior_n + o.data_n) * (mean_diff'*mean_diff);
            [~, p] = cholcov(o.post_cov);
            assert(p == 0, 'o.posterior covariance was not PD!');
            
            % Murphy (228)
            o.pred_dof = o.post_dof - o.dim + 1;
            o.pred_cov = (o.post_n + 1) / (o.post_n * o.pred_dof) * o.post_cov;
            o.pred_mvtparams = make_mvt(o.post_mean, o.pred_cov, o.pred_dof);
        end                       
        
        function [p] = pred_like(o, X)
% like(X) returns the posterior predictive likelihood of X
% TODO: Figure out the API for the non-conjugate case. Store the state as
% fields?
%
% Really need a distribution class.
            p = fast_mvtpdf(o.pred_mvtparams, X);
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
