function [ like ] = beta_intensity_like( mm, k, x )
% Beta intensity likelihood for independent mixture model.
%
% mm.beta_intensity_prior - [d r k]; see beta_prior_updf.
% mm.s_z                  - class indicator
% mm.K                    - number of classes
% k                       - class we are computing for
% x                       - column vector of intensities

    [Nk, D] = size(x);
    assert(D == 1, 'x must be a column');
    % Parameters
    Nsamps = 50;
    burn_in = 100;
    MAX_LIKE = 1e10;
    
    % Update step: compute the posterior     
    like = zeros(Nk, 1);
    [mm.beta_posterior_samps{k} accept_rate] = ...
        beta_posterior_samples( mm.prior_beta_params, x, Nsamps, burn_in, 0.01 );    
    
    % Data-likelihood step.
    % Compute a kernel density estimate to mm.beta_posterior_samps{k}.
    % Then, evaluate the data likelihood for beta distributions picked at
    % evenly spaced points on this kernel density.
    %
    % Actually, make this a TODO. Just evaluate each sample. We should get
    % the probabilistic behavior by drawing many, many distribution, and
    % the general smoothness of the posterior.
    s = mm.beta_posterior_samps{k}(:,1);                          
    m = mm.beta_posterior_samps{k}(:,2);
    [a b] = beta_sm_to_ab(s, m);        
    
    for ns = 1:Nsamps
        like = like + betapdf(x, a(ns), b(ns));
    end

    like = like / Nsamps;
    
    disp(['Accept rate: ', num2str(accept_rate)]);
    figure(10);
    subplot(2,1,1);    
    hist(s);
    title('s samples');
    
    subplot(2,1,2);
    hist(m);
    title('m samples');
    drawnow;
    
    % Clamp the top
    like(like == Inf) = MAX_LIKE;
    
    assert(all(like >= 0 & like ~= Inf));    
end

