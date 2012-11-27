function [ samples ] = run_mdp( iter_f, state, params, X, img, mdp, niters, burn_in, sample_freq )

    loglike = zeros(niters, 1);    
    img_h = figure;
    set(img_h, 'Units', 'normalized', 'position', [0.1 0.1 0.9 0.4]);
    
    for n = 1:niters
        tic
        
        cluster_assigns_old = mdp.cluster_assigns;
        [mdp, state, ll] = iter_f(mdp, state, params);
        
        % Timing and diagnostics
        change = mean(abs(cluster_assigns_old ~= mdp.cluster_assigns));        
        itertime = toc;        
        rate     = mdp.N / itertime;
        
        fprintf(1, 'Iter %d: time = %0.2f; speed = %.2f; K = %d, change = %0.3f; loglike = %06.2f\n', ...
                n, itertime, rate, mdp.n_clusters, change, ll);

        % Save the log-likelihood
        loglike(n) = ll;                
        
        if mod(n, sample_freq) == 0
            plot_img(img_h, n, loglike, X, img, mdp);
            %plot_mdp_diag(diag_h, mdp, X);
           
            drawnow;
            
            if n > burn_in
                % TODO:
            end            
        end
    end
end

