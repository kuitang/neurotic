function [ samps, accept_rate ] = beta_posterior_samples( prior_params, x, Nsamps, Nburn_in, stdev )
% [s m] = Nsamps x 2 matrix of variance-multipliers and means from beta
% posterior.
%
% [d r k] = prior_params (see beta_prior_updf)
% x       = vector of observations
% Nsamps  = number of samples desired
% burn_in = number of initial samples to discard
% stdev   = standard deviation for lognormal proposal for M-H (try 0.1 or
%           0.01)
%
    % Sample!
    samps = zeros(Nsamps, 2);
    
    % MLE init is stupid and underestimates mean. Do method of moments
    % instead.
    m = mean(x);
    m_star = interval_to_positive(m);
    % Notebook pp 29
    s = var(x) / (m * (1 - m));
            
    % Our prior is designed for s > 2.
    % Push away for a low-likelihood region
    s = max(4, s);
        
    for n = 1:Nburn_in
        [s m_star] = beta_posterior_mh(prior_params, s, m_star, x, stdev);
    end

    samps(1,:) = [s m_star];

    Naccept = 0;    
    for n = 2:Nsamps
        s_old      = samps(n - 1,1);
        m_star_old = samps(n - 1,2);
        [s_new m_star_new acc] = beta_posterior_mh(prior_params, s_old, m_star_old, x, stdev);
        samps(n,:) = [s_new m_star_new];
        if acc        
            Naccept = Naccept + 1;
        end
    end
    
    accept_rate = Naccept / Nsamps;
end
