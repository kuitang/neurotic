function [ up ] = beta_prior_updf( s, m, d, r, k )
% Prior for mean/variance-scale parameterization of Beta.
% Equation Bouguila (18), Roberts (7).
%
% p - unnormalized probability density
% s := a + b; m := a / s where a, b of standard parameterization
%
% d - (2, 0.5) distance penalty
% r - skew penalty
% k - variance penalty
%
% Bouguila uses [d r k] = [3 0.1 0.00001]

    log_t1 = log((1 - exp(-d * ((s - 2).^2 + (m - 0.5).^2))));
    log_t2 = -r / (s.^2 .* m .* (1 - m)) - (k * s.^2) / 2;
    
    up = exp(log_t1 + log_t2);
    assert(all(~isinf(up) & ~isnan(up) & up >= 0));        

end

