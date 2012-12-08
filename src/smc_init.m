function [ particles, weights ] = smc_init( mdp, background_like )
% Runs one SMC pass through X to obtain sensible initialization for
% subsequent MCMC steps.
%
% mdp - an uninitialized MDP object (mdp.cluster_assigns = zeros(...))
%       Returns a fully populated mdp.cluster_assigns as well as
%       mdp.cluster_likes.
% X   - N x D data matrix
%
% References
% [1] Paul Fearnhead, Particle filters for mixture models with an unknown
%     number of components, Statistics and Computing 14: 11-21, 2004.
% [2] Paul Fearnhead and Peter Clifford, On-line inference for hidden
%     Markov models via particle filters, J. R. Statist. Soc. B (2003) pp.
%     887-899
% [3] Frank Wood and Michael Black, A nonparametric Bayesian alternative
%     to spike sorting, Journal of Neuroscience Methods, 2008,
%     doi:10.1016/j.jneumeth.2008.04.030,
    
    X = mdp.X;
    [N, D] = size(X);           
    
    % Preallocate
    
    % Maintain a sub particle vector and super particle matrix with
    % weights.
    Nsub = 20;
    MAX_CLUSTERS = 100;
    
    particles(1, Nsub) = MDP;
    weights(1, Nsub) = 0;    
        
    % Intialize
    %
   
    perm = randperm(N);
    
    % IMPORTANT! We have to start off by giving mdp a background class.    
    mdp.cluster_likes{1} = background_like;
    % Now MANUALLY initialize. We can't use mdp.assign, because that would
    % create a new "normal" cluster.
    mdp.n_clusters = 1;
    mdp.cluster_assigns = zeros(N, 1);
    mdp.cluster_assigns(1) = 1;
    mdp.cluster_counts(1) = 1;
    
    for p = 1:Nsub
        particles(p) = copy(mdp);
    end   
    weights(:)   = 1 / Nsub;    

    % Precompute, as in neal3
    prior = copy(mdp.prior);
    prior.fit([]);
    prior_pred_like = prior.pred_like(X); % col vector
    
    super_shape = [MAX_CLUSTERS + 1, Nsub];
    tic
    for nn = 2:N
        
        if mod(nn, 1000) == 0
            ktime = toc;
            speed = 1000 / ktime;
            disp(['Iter ' num2str(nn) ' time/k = ' num2str(ktime) ' speed = ' ...
                  num2str(speed) ' /s']);
            tic;
        end
        
        n = perm(nn);   
        % Clear junk from all previous iterations
        super_particles(super_shape(1), super_shape(2)) = MDP;
        super_weights = zeros(super_shape);        
        
        % EVOLUTION STAGE
        % For each particle, evolve to each possible next state.
        % Fearnhead (2003) pp. 890        
        for p = 1:Nsub
            % Generate a transition probability table.
            % transition(r, c) = P(class r | class c)
            pmdp = particles(p);
            transition = [ pmdp.cluster_counts ; pmdp.concentration ] ./ ...
                         ( pmdp.N - 1 + pmdp.concentration );
            
            assert(pmdp.n_clusters <= MAX_CLUSTERS, 'MAX_CLUSTERS set too low')
            for k = 1:(pmdp.n_clusters + 1)
                super_particles(k,p) = copy(particles(p));
                
                super_particles(k,p).assign(n, k);
                super_particles(k,p).N = super_particles(k,p).N + 1;                               
                
                if k == pmdp.n_clusters + 1
                    % New class
                    like = super_particles(k,p).cluster_likes{k}.pred_like_scalar(X(n,:));
                else
                    like = prior_pred_like(n);
                end
                
                super_weights(k,p) = weights(p) * transition(k) * like;                                                     
            end    
        end
        
        % RESAMPLING STAGE
        [new_flat_sub_weights, c] = fearnhead_resample(super_weights(:), Nsub);
        new_sub_weights = reshape(new_flat_sub_weights, super_shape);
        new_nnzidxs = new_sub_weights > 0;
        n_nnz = sum(new_nnzidxs(:));
        
        % n_nnz may be off by Nsub by 1.
        assert(abs(n_nnz - Nsub) <= 1, 'fearnhead_resample did not return correct sized sample');
        if n_nnz < Nsub
            % Add in a particle that was zero
            [i, j] = find(~new_nnzidxs, 1);            
            new_sub_weights(i, j) = 1/c;
            new_nnzidxs(i, j) = 1;
        elseif n_nnz > Nsub
            % Remove a particle that was 1/c
            [i, j] = find(new_sub_weights == 1/c, 1);
            new_sub_weights(i, j) = 0;
            new_nnzidxs(i, j) = 0;
        end
        
        % Subindex
        particles = reshape(super_particles(new_nnzidxs), 1, Nsub);
        weights   = reshape(new_sub_weights(new_nnzidxs), 1, Nsub);
    end
end
