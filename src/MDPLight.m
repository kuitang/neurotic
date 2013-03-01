classdef MDPLight < matlab.mixin.Copyable
%MDP Only Bookkeeping for Mixture of Dirichlet Process distribution.
% This version doesn't care about likelihoods.
%
% Also supports adding arbitrarily many new classes.
%
% References
% [1] Radford Neal, Markov Chain Sampling Methods for Dirichlet Process
%     Mixture Models
    
    properties
        concentration
        
        % Just storage
        N = 0;
        data_dim
        X
        
        misc_data        
        % prior is a struct with fields specified by the particular user of
        % this class.
        prior
        
        % This is part of the Markov chain (as are the cluster
        % assignments), but we need to manage the sampled parameters (put
        % in state) together with the rest of the clusters for
        % adding/removing/etc.
        state
        
        % N x 1 vector with range [0, num_clusters]. 0 represents
        % unassigned and 1 represents background.
        cluster_assigns = [];
        % Number of currently active clusters
        n_clusters = 0;
        % Summary statistics
        cluster_counts = [];                        
    end
    
    methods
        function o = MDPLight(concentration_or_rhs, X, prior, state, misc_data, cluster_assigns)
            if nargin == 1 % copy constructor
                % Assign the arrays
                rhs = concentration_or_rhs;
                o.concentration = rhs.concentration;
                o.N = rhs.N;
                o.data_dim = rhs.data_dim;
                o.X = rhs.X;
                o.misc_data = rhs.misc_data;
                o.prior = rhs.prior;
                o.state = rhs.state;
                
                o.cluster_assigns = rhs.cluster_assigns;
                o.n_clusters = rhs.n_clusters;
                o.cluster_counts = rhs.cluster_counts;
                
                % INITIALIZE THE STATE
                
            elseif nargin > 0 % support empty constructor for array construction
                % cluster_assigns must be all specified, for now.            
                [o.N, o.data_dim] = size(X);
                o.X = X;
                o.misc_data = misc_data;
                o.concentration = concentration_or_rhs;
                o.misc_data = misc_data;                
                o.prior = prior;     
                o.state = state;
                
                nzidxs = cluster_assigns > 0;
                if sum(nzidxs > 0)                
                    % Initialize our model with cluster_assigns
                    o.cluster_assigns = cluster_assigns;            
                    o.n_clusters = max(cluster_assigns);                    
                    o.cluster_counts = accumarray(cluster_assigns(nzidxs), ones(o.N, 1));
                end
            end
        end
                        
        function k_old = remove_point(o, n)
% k_old = remove_point(o, n) removes n. First step in Gibbs sampling.
            k_old = o.cluster_assigns(n);
            assert(k_old ~= 0, 'Point was never assigned in the first place');
            o.cluster_counts(k_old) = o.cluster_counts(k_old) - 1;
            o.cluster_assigns(n) = 0;                                    
        end
        
        function assign(o, n, k)
% assign(o, n, k) assigns n to k in the general case.
% After all assignments are done, you must call kill_empty_clusters!            
            if k > o.n_clusters
                % New cluster                
                % Implicitly resize the count array
                o.cluster_counts(k)  = 1;                
                o.cluster_assigns(n) = k;   
                o.n_clusters = k; % But not really, if we inserted gaps
            else
                % Existing                
                o.cluster_counts(k) = o.cluster_counts(k) + 1;
                o.cluster_assigns(n) = k;
            end                            
        end
                
        function kill_cluster(o, k)            
            % To delete an element in a cell array, use () indexing.
            o.cluster_likes(k) = [];
            o.cluster_counts(k) = [];
            o.state(k) = [];
            
            % For now, we must have an empty class to begin with
            assert(sum(o.cluster_assigns == k) == 0, 'unassigning unsupported');            
            
            % Collapse
            gtidxs = o.cluster_assigns > k;
            o.cluster_assigns(gtidxs) = o.cluster_assigns(gtidxs) - 1;
            o.n_clusters = o.n_clusters - 1;
            
            % Update the my_k's
            for kk = k:o.n_clusters
                if kk > 1
                    assert(o.cluster_likes{kk}.pdfs{2}.my_k - 1 == kk, ...
                           'GraphDist my_k precondition failed');
                    o.cluster_likes{kk}.pdfs{2}.my_k = kk;
                end
            end
            
            assert(length(o.cluster_counts) == o.n_clusters && ...
                   length(o.cluster_likes) == o.n_clusters, 'bad removal in kill_cluster');
        end
        
        % DESTROY DEAD CLUSTERS!
        function prune_clusters(o, k)
            for k = 1:o.n_clusters
                if o.cluster_counts(k) == 0
                    assert(k ~= 1, 'uh oh, killing the background class!');
                    o.kill_cluster(k);
                end
            end
            
            % In addition, kill state that we don't have
            extra_states = (o.n_clusters + 1):length(o.state);
            o.state(extra_states) = [];
        end                        
    end
end

