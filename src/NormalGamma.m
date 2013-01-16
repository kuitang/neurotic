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
        pred_precision;                
    end
    
    properties (Dependent)
        pred_mean, pred_dof;
    end            
    
    methods
        function pm = get.pred_mean(o)
            pm = o.post_mean;
        end                
        
        function pd = get.pred_dof(o)
            % Murphy (110)
            pd = 2 * o.post_shape;
        end
        
        function o = NormalGamma(prior_mean_or_rhs, prior_n, prior_shape, prior_rate)
            if nargin == 1 % copy constructor
                % REMEMBER ME!!! Standard trick to copy fields. 
                fns = properties(prior_mean_or_rhs);
                for i=1:length(fns)
                    o.(fns{i}) = rhs.(fns{i});
                end
            else    
                o.prior_mean  = prior_mean_or_rhs;
                o.prior_n     = prior_n;
                o.prior_shape = prior_shape;
                o.prior_rate  = prior_rate;                                
                
                % IMPORTANT: Make the distribution ready for prediction and
                % update
                o.fit([]);
            end
        end                
                
        function fit(o, x)            
% fit(X) computes posterior and predictive parameters. x vector.
            [o.data_n, dim] = size(x);
            
            % Special case for empty data
            if o.data_n == 0                
                data_mean = 0;
                data_cov  = 0;
            elseif o.data_n == 1
                data_mean = x;
                data_cov  = 0;                            
            else 
                assert(dim == 1);            
                data_mean = mean(x);
                data_cov  = cov(x, 1);            
            end            
                                               
            % Murphy (85) to (89)
            o.post_n     = o.prior_n + o.data_n;
            o.post_shape = o.prior_shape + o.data_n / 2;
            
            o.post_mean  = (o.prior_n*o.prior_mean + o.data_n*data_mean) / ...
                           (o.prior_n + o.data_n);
            
            % Rate is sum of squares                        
            
            o.post_rate  = o.prior_rate + o.data_n / 2 * data_cov + ...
                           o.prior_n*o.data_n*(data_mean - o.prior_mean)^2 ./ ...
                           (2*(o.prior_n + o.data_n));
            
            % Murphy (110)
            o.pred_precision = (o.post_shape * o.post_n) / (o.post_rate * (o.post_n + 1));                       
        end        
        
        % Online updates. Currently only online updates for mean and cov.
        % Eventually need to update mvtparams...
        function add_point(o, x)
            % CURRENTLY NOT SUPPORTED
        end
        
        function remove_point(o, x)
            % CURRENTLY NOT SUPPORTED
        end
        
        function [p] = pred_like(o, x)
% like(X) returns the posterior predictive likelihood of X
            % Standardize
            z = (x - o.post_mean) .* sqrt(o.pred_precision);
            p = tpdf(z, o.pred_dof);
        end
        
        function [p] = pred_like_scalar(o, x)
            p = o.pred_like(x);            
        end
        
        function [x] = sample_prior_pred(o, n)
% x = sample_prior(n) samples one mean and precision and n iid datapoints
            precision = gamrnd(o.prior_shape, 1 / o.prior_rate);
            prior_sd  = sqrt((o.prior_n * precision)^-1);
            mu        = normrnd(o.prior_mean, prior_sd);
            
            % Murphy (72)
            z = trnd(2 * o.prior_shape, n, 1);
            v = o.prior_rate ./ (o.prior_shape * o.prior_n);
            x = sqrt(v) .* z + mu;
        end
        
        function [x] = sample_posterior_pred(o, n)
% X = sample_posterior(n) samples one mean and cov and n iid datapoints
            precision = gamrnd(o.post_shape, 1 / o.post_rate);
            post_sd   = sqrt((o.post_n * precision)^-1);
            mu        = normrnd(o.post_mean, post_sd);
            
            % Murphy (72)
            z = trnd(2 * o.post_shape, n, 1);
            v = o.post_rate ./ (o.post_shape * o.post_n);
            x = sqrt(v) .* z + mu;            
        end
    end
    
end
