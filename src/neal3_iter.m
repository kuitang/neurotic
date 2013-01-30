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
    
    % Hack: We haven't assigned my_k, so we can't use this
    prior_pred_like = prior.pdfs{1}.pred_like(mdp.X(:,1:4));
    % A new class has distance 0 from its cluster center.
    prior_pred_like = prior_pred_like .* prior.pdfs{2}.underlying.pred_like(0);
    
    x_max = max(mdp.X(:,1));
    y_max = max(mdp.X(:,2));
    
    ll = 0;
    
    actual_n = 0;
    
    for n = randperm(mdp.N)        
        k_old = mdp.remove_point(n);                        
        
        z_pdf = [ mdp.cluster_counts ; mdp.concentration ] ./ ...
                    ( mdp.N - 1 + mdp.concentration );
        x_pdf  = zeros(size(z_pdf));
        
        for k = 1:mdp.n_clusters
            x_pdf(k) = mdp.cluster_likes{k}.pred_like_scalar(mdp.X(n,:));
        end
        
        % New cluster
        x_pdf(end) = prior_pred_like(n);                
        
        % HACK
%         for k = 1:mdp.n_clusters
%             X5_before = mdp.X(:,5);
%             mdp = refit_graph_dist(mdp, x_max, y_max, k);
%             assert(length(mdp.misc_data.slic_inds{end}) > 1);
%             x_pdf(k) = mdp.cluster_likes{k}.pred_like_scalar(mdp.X(n,:));            
%             
%             if k > 1 && mod(actual_n, 1000) == 0
%                 h = figure(30);
%                 sh = subplot(2, mdp.n_clusters, 1);
%                 imshow(showslic(mdp.misc_data.segments, mdp.misc_data.img));
%                 
%                 plot_graph_edge_dist(h, k, mdp, x_max, y_max);
%             end
%         end        
        
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

