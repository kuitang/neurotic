function [ up ] = beta_prior_updf( params, s, m, logp )
% Prior for mean/variance-scale parameterization of Beta.
% Equation Bouguila (18), Roberts (7).
%
% Modified to force s > 2. We want a unimodal beta, not one with diverging
% at 0 or 1.
%
% up        = unnormalized probability density
% s         = a + b (variance multiplier)
% m         = a / s (mean)
%
% [ d r k ] = scalar params where
%
% d         = distance to (2, 0.5) preference
% r         = skew penalty
% k         = variance penalty
%
% log       = [optional] output log probability
%
% Bouguila uses [d r k] = [3 0.1 0.00001]

    d = params(1); r = params(2); k = params(3);
    
%     % Box-constrained to a, b > 1.
     [a b] = beta_sm_to_ab(s, m);
%     if a < 1 || b < 1    
%         up = -Inf;
%     else
        log_t0 = log(1 - exp(-( (s - 0.5).^2))) + log(1 - exp(-( (a - 1).^2))) + log(1 - exp(-( (b - 1).^2)));
        log_t1 = log((1 - exp(-d * ((s - 2).^2 + (m - 0.5).^2))));
        log_t2 = -r / (s.^2 .* m .* (1 - m)) - (k * s.^2) / 2;

        up = log_t0 + log_t1 + log_t2;   
%     end
    if nargin > 3 && logp
        return
    end
    
    up = exp(up);
    assert(all(~isinf(up) & ~isnan(up) & up >= 0));        

end

