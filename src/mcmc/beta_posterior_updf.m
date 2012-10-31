function [ up ] = beta_posterior_updf( prior_params, s, m, x, logp )
% Posterior of the beta_prior_updf distribution.
% Bouguila (19) and (20)
%
% Used as a subroutine in Beta mixture models. Class labels don't exist
% here: handle them (through indexing the big X) before calling this.
%
% up      = unnormalized probability density
%
% [d r k] = prior_params (see beta_prior_updf)
% s, m    = variance multiplier and mean
% K       = number of classes
% logp    = [optiona;] output log probability

    [N D] = size(x);
    assert(D == 1, 'beta distribution is for scalar x!');
    assert(length(s) == 1 && length(m) == 1, 's and m must be scalar!');        

    log_prior = beta_prior_updf(prior_params, s, m, true);
    
    sm  = s * m;
    s1m = s * (1 - m);
    log_Z     = N * (gammaln(s) - gammaln(sm) - gammaln(s1m));        
    log_like  = sm * sum(log(x)) + s1m * sum(log1p(-x));
    up = log_prior + log_Z + log_like;
    
    if nargin > 4 && logp
        return
    end
    
    up = exp(up);
    assert(all(~isinf(up) & ~isnan(up) & up >= 0));

end

