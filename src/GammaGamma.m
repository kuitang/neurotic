classdef GammaGamma < OnlineDistribution
% Gamma likelihood with conjugate Gamma **rate** prior.
%
% This parameterization:
%
%   p(x | shape, rate) = rate^shape / gamma(shape) * x^(shape - 1) * 
%                        exp(-rate * x)
%
% The posterior update equations are well known (i.e. Wikipedia). The
% predictive likelihood compound gamma, which simplifies to beta. See [1].
%
% References
% [1] Satya D. Dubey, Compound Gamma, Beta, and F Distributions, Metrika
%     1970: 27-31.
    
    properties
        dims;
        shape, prior_shape, prior_rate;
        post_shape, post_rate;
        data_n;        
    end
    
    methods
        function o = GammaGamma(shape, prior_shape, prior_rate)
            % dims - the dims of the data to capture
            if nargin == 1 % copy constructor
                % REMEMBER ME!!! Standard trick to copy fields. 
                fns = properties(prior_mean_or_rhs);
                for i=1:length(fns)
                    o.(fns{i}) = rhs.(fns{i});
                end                
            end            
            o.data_n = 0;
            o.shape = shape;
            o.prior_shape = prior_shape;
            o.prior_rate  = prior_rate;
            
            o.post_shape  = prior_shape;
            o.post_rate   = prior_rate;
        end
        
        function fit(o, X)
            o.data_n = length(X);
            o.post_shape = o.prior_shape + o.data_n * o.shape;
            o.post_rate  = o.prior_rate  + sum(X);
            assert(o.post_rate > 0);
        end
        
        function add_point(o, x)
            o.data_n = o.data_n + 1;
            o.post_shape = o.post_shape + o.shape;
            o.post_rate  = o.post_rate  + x;
            assert(o.post_rate > 0);
        end
        
        function remove_point(o, x)
            % If we're going to zero data, revert to the prior
            if o.data_n <= 1
                o.fit([]);
            else            
                o.data_n = o.data_n - 1;
                o.post_shape = o.post_shape - o.shape;
                o.post_rate  = o.post_rate - x;
                assert(o.post_rate > 0);
            end            
        end
        

        
        % Double check the paper with this
        function xs = sample_prior(o, n)
            rate = gamrnd(o.prior_shape, 1 / o.prior_rate);
            xs = gamrnd(o.shape, 1 / rate, n, 1);
        end
        
        function p = pred_like_scalar(o, X)            
            % Dubey (1.1)
            if X == 0 % avoid log(0)
                log_top_part = o.post_shape * log(o.post_rate);
                log_bot_part = betaln(o.shape, o.post_shape) + (o.shape + o.post_shape) * log(o.post_rate + X);
                p = exp(log_top_part - log_bot_part);
                % Important
                p = p * 0^(o.shape - 1);
            else
                log_top = o.post_shape * log(o.post_rate) + (o.shape - 1) * log(X);
                log_bot = betaln(o.shape, o.post_shape) + (o.shape + o.post_shape) * log(o.post_rate + X);
            
                p = exp(log_top - log_bot);                  
            end
        end
        
        function p = pred_like(o, X)
        % Predictive likelihood (compound gamma). X can be a vector.
            % Dubey (1.1)
            p = o.pred_like_scalar(X);            
        end
    end
    
end

