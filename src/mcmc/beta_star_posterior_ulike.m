function [ up ] = beta_star_posterior_ulike( prior_params, s, m_star, x )
% Log likelihood of the posterior.
% Like beta_posterior_updf, with m_star = m / (1 - m)

% Bouguila bottom paragraph on page 220
    
    %e_m_star = exp(m_star);
    
    t = positive_to_interval(m_star);
    log_deriv = -2 * log(1 + m_star);    
            
    up = beta_posterior_updf(prior_params, s, t, x, true) + log_deriv;

end

