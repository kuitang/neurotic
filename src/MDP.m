classdef MDP < handle
%MDP Bookkeeping for Mixture of Dirichlet Process distribution. Use in
% Gibbs algorithms.
%
% Includes special background class #1.
%
% References
% [1] Radford Neal, Markov Chain Sampling Methods for Dirichlet Process
%     Mixture Models
    
    properties
        concentration
        
        % Just storage
        N
        data_dim
        X
        prior
        
        % N x 1 vector with range [0, num_clusters]. 0 represents
        % unassigned and 1 represents background.
        cluster_assigns
        % Number of currently active clusters
        n_clusters
        % Summary statistics
        cluster_counts
        
        background_like
        cluster_likes
    end
    
    methods
        function o = MDP(concentration, X, cluster_assigns, background_like, prior)
% cluster_assigns must be all specified, for now.            
            [o.N, o.data_dim] = size(X);
            o.X = X;
            o.concentration = concentration;            
            o.background_like = background_like;        
            o.prior = prior;
            
            % Initialize our model with cluster_assigns
            o.cluster_assigns = cluster_assigns;            
            nzidxs = cluster_assigns > 0;
            o.n_clusters = max(cluster_assigns);
            assert(sum(nzidxs) == o.N, 'uninitialized start not implemented');
            o.cluster_counts = accumarray(cluster_assigns(nzidxs), ones(o.N, 1));            
            
            % TODO: Perhaps remove the restriction
            
                        
            if o.n_clusters <= 1
                % Always create a "dummy" cluster entry to account for
                % background
                o.cluster_likes = cell(1, 1);
            else
                o.cluster_likes = cell(o.n_clusters, 1);
            end
            
            if isa(background_like, 'matlab.mixin.Copyable')
                o.cluster_likes{1} = copy(background_like);
            else
                o.cluster_likes{1} = background_like;
            end
            
            for k = 2:o.n_clusters
                idxs = cluster_assigns == k;
                
                o.cluster_likes{k} = copy(prior);
                o.cluster_likes{k}.fit(X(idxs,:));
            end
        end
        
        function [out] = light_copy(o)
            error('not implemented!');
        end
        
        function refit(o, k)
% Refits the likelihood for k for all points currently assigned to k
%
% TODO: Online updates (should integrate with remove_point and
% assign_to_new

            % If we are not background...
            if k > 1
                idxs = o.cluster_assigns == k;
                o.cluster_likes{k}.fit(o.X(idxs,:));
            end
        end
        
        function k_old = remove_point(o, n)
% Remove item n. First step in Gibbs sampling. Refits the distribution.
            k_old = o.cluster_assigns(n);
            o.cluster_counts(k_old) = o.cluster_counts(k_old) - 1;
            o.cluster_assigns(n) = 0;
            
            % TODO: Online updates.
            if k_old > 1
                o.cluster_likes{k_old}.remove_point(o.X(n,:));
            end
            %o.refit(k_old);
        end
        
        function assign(o, n, k)
% assign(o, n, k) assigns n to k in the general case. If k does not exist,
% it is created with a new prior (now posterior).
            assert(k < o.n_clusters + 2, 'k is too high!');
            if k > o.n_clusters
                % New
                dist = copy(o.prior);                
                o.assign_to_new(n, dist);
                dist.fit(o.X(n,:));
            else
                % Existing
                o.assign_to_existing(n, k);
                if k > 1
                    o.cluster_likes{k}.add_point(o.X(n,:));
                end
                %o.refit(k);
            end                            
        end
        
        function assign_to_new(o, n, cluster_like)
% add_cluster(n, cluster_like) assigns point n to a new cluster whose
% distribution is cluster_like.
% If you're Gibbs sampling, make sure to call REMOVE first!            
    
            o.n_clusters = o.n_clusters + 1;
            o.cluster_counts = [ o.cluster_counts ; 1 ];
            % Need to do this to prevent collapsindbg heterogenous members
            % (otherwise, the background class disappears and you only get
            % the NormalWisharts)
            o.cluster_likes =  [ o.cluster_likes ; 1 ];
            o.cluster_likes{end} = cluster_like;            
            o.cluster_assigns(n) = o.n_clusters;
            
            assert(length(o.cluster_counts) == o.n_clusters && ...
                   length(o.cluster_likes)  == o.n_clusters, ...
                   'bad augmentation in add_cluster');
        end

        function assign_to_existing(o, n, k)
% This DOES NOT update posteriors! You have to do that yourself.
% If you're Gibbs sampling, make sure to call REMOVE first!

            assert(k > 0 && k <= o.n_clusters, 'k must exist!');
                        
            o.cluster_counts(k) = o.cluster_counts(k) + 1;
            o.cluster_assigns(n) = k;
        end    
        
        function like = kill_cluster(o, k)
% like = kill_cluster(k) places cluster k into like and removes from our
% list. You are responsible for freeing memory if you really want to get
% rid of it.            
            like = o.cluster_likes;
            
            % To delete an element in a cell array, use () indexing.
            o.cluster_likes(k) = [];
            o.cluster_counts(k) = [];
            
            % For now, we must have an empty class to begin with
            assert(sum(o.cluster_assigns == k) == 0, 'unassigning unsupported');            
            
            % Collapse
            gtidxs = o.cluster_assigns > k;
            o.cluster_assigns(gtidxs) = o.cluster_assigns(gtidxs) - 1;
            o.n_clusters = o.n_clusters - 1;
            
            assert(length(o.cluster_counts) == o.n_clusters && ...
                   length(o.cluster_likes) == o.n_clusters, 'bad removal in kill_cluster');
        end
                
    end
end

