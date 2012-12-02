function [ new_weights, c ] = fearnhead_resample( weights, n_particles )
% new_weights = fearnhead_resample( weights, n_particles )
% Fearnhead's exhaustive resampling N < RN algorithm.
%
% weights     - M-vector of weights (full weights)
% n_particles - < M to resample
%
% new_weights - M-vector of weights, only n_particles nonzero.
%
% References
% [1] Paul Fearnhead and Peter Clifford, On-line inference for hidden
%     Markov models via particle filters, J. R. Statist. Soc. B (2003) pp.
%     887-899
    
    M = length(weights);
    Mones = ones(size(weights));

    assert(n_particles < M, 'frs:tooManyParticles', 'n_particles > length(weights); problem.');
    assert(n_particles > 0, 'frs:zeroParticles', 'n_particles must be > 0');
            
    % TODO: Analytic form? \sum min(cq_j, 1) = N
    obj = @(c) sum(min(c * weights, Mones)) - n_particles;    
    c = fsolve(obj, 1, optimset('display', 'off'));    
    
    bigidxs = weights > 1/c;
    n_big = sum(bigidxs);
    
    n_small = n_particles - n_big;
    assert(n_small > 0, 'frs:noSmall', 'constraint n_particles - n_big > 0 violated!');
        
    smallidxs = stratified_resample(weights, ~bigidxs, n_small);
    
    new_weights = zeros(size(weights));
    new_weights(bigidxs)  = weights(bigidxs);
    new_weights(smallidxs) = 1/c;    
end


