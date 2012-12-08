function [ mdp, state, ll ] = neal8_iter( X, mdp, state, params )
% Gibbs sampling iteration for Algorithm 8 in Neal 2000 [1]
%
% state: 
% References:
% [1] Radford Neal, Markov Chain Sampling Methods for Dirichlet Process
%     Mixture Models
    
    unrealized_likes = cell(params.m, 1);
    
    ll = 0;
    for n = 1:mdp.N
        k_old = mdp.remove_point(n);
        
        if mdp.cluster_count(k_old) == 0
            % Singleton: consider the current distribution as a new
            % candidate
            m = params.m - 1;
            unrealized_likes(params.m) = mdp.cluster_likes{k_old};
        else
            m = params.m;
        end
        
        % Draw m samples from the prior
        
        
    end    

end

