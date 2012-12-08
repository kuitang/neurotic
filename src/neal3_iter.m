function [ mdp, state, ll ] = neal3_iter( mdp, state, params )
% Gibbs sampling iteration for Algorithm 3 in Neal 2000 [1]
%
% state: 
% References:
% [1] Radford Neal, Markov Chain Sampling Methods for Dirichlet Process
%     Mixture Models

% Precompute
    prior = copy(mdp.prior);
    prior.fit([]);
    prior_pred_like = prior.pred_like(mdp.X); % col vector
    
    x_max = max(mdp.X(:,1));
    y_max = max(mdp.X(:,2));
    
    ll = 0;
    
    actual_n = 0;
    
    for n = randperm(mdp.N)

        
        k_old = mdp.remove_point(n);        
        
        % Removed a point, so we must refit the posterior
        %
        % If our model is fully online, not a problem.
        %mdp.refit(k_old);
        
        z_pdf = [ mdp.cluster_counts ; mdp.concentration ] ./ ...
                    ( mdp.N - 1 + mdp.concentration );
        x_pdf  = zeros(size(z_pdf));        
        % New cluster
        x_pdf(end) = prior_pred_like(n);                
        
        for k = 1:mdp.n_clusters
            refit_graph_dist(mdp, x_max, y_max, k_old);
            x_pdf(k) = mdp.cluster_likes{k}.pred_like_scalar(mdp.X(n,:));
        end
        
        pdf = z_pdf .* x_pdf;
        
        % Unnormalized inverse cdf sampling        
        cdf = cumsum(pdf);        
        k_new = find(cdf > cdf(end)*rand(1), 1);
        mdp.assign(n, k_new);        
        
        % DESTROY DEAD CLUSTERS!
        if mdp.cluster_counts(k_old) == 0
            assert(k_old ~= 1, 'uh oh, killing the background class!');
            mdp.kill_cluster(k_old);
        end
        
        actual_n = actual_n + 1;
        if mod(actual_n, 1000) == 0
            disp(['actual_n = ' num2str(actual_n)]);
            disp(['x_pdf = ' num2str(x_pdf')]);
        end
        
        ll = ll + log(pdf(k_new));
    end

end

