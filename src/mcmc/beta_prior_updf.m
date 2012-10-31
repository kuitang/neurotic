function [ up ] = beta_prior_updf( params, s, m, logp )
% Prior for mean/variance-scale parameterization of Beta.
% Equation Bouguila (18), Roberts (7).
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
    
    log_t1 = log((1 - exp(-d * ((s - 2).^2 + (m - 0.5).^2))));
    log_t2 = -r / (s.^2 .* m .* (1 - m)) - (k * s.^2) / 2;
    
    up = log_t1 + log_t2;    
    if nargin > 3 && logp
        return
    end
    
    up = exp(up);
    assert(all(~isinf(up) & ~isnan(up) & up >= 0));        

end

