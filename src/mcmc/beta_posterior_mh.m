function [ s, m_star, accept ] = beta_posterior_mh( prior_params, s_old, m_star_old, x, stdev )
% Sample from the beta posterior. Metropolis-Hasting step.
% Bouguila pages 219 and 220
%
% s_new, m_new = next (or recycled) samples from the Markov chain.
%
% [d r k]      = prior_params (see beta_prior_pdf)
% s_old, m_old = previous samples from the Markov chain
% x            = data
% stdev        = standard deviation of log-normal proposal
%                (Bouguila uses 0.1 = sqrt(0.01))
            
    % Draw samples
    % Bouguila page 220
    s      = lognrnd(log(s_old),      stdev);
    m_star = lognrnd(log(m_star_old), stdev);
    accept = true;
        
    % Compute M-H acceptance criterion
    % Note MATLAB's lognlike returns NEGATIVE log-likelihood, but we want
    % the positive.
    log_top = beta_star_posterior_ulike(prior_params, s, m_star, x) - ...
              lognlike([log(m_star) stdev], m_star_old) - ...
              lognlike([log(s)      stdev], s_old);
      
    log_bot = beta_star_posterior_ulike(prior_params, s_old, m_star_old, x) - ...
              lognlike([log(m_star_old) stdev], m_star) - ...
              lognlike([log(s_old)      stdev], s);                  
    
    r = exp(log_top - log_bot);
    assert(~isnan(r));    
    
    if rand > r % M-H reject        
        s = s_old;
        m_star = m_star_old;
        accept = false;
    end        

end
