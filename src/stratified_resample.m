function [ sampleidxs ] = stratified_resample( weights, subsetidxs, n_particles )
% Resample n_particles from weights(subset). Actual samples may differ from
% n_particles by at most 1 [2].
%
% weights     - N-vector of all weights
% subsetidxs  - N-vector indexing the set to resample from
% n_particles - number of particles to select
% sampleidxs  - N-vector (n_particles nonzero) indexing the resample
%
% References
% [1] Paul Fearnhead and Peter Clifford, On-line inference for hidden
%     Markov models via particle filters, J. R. Statist. Soc. B (2003) pp.
%     887-899
% [2] J. Carpenter, P. Clifford, and P. Fearnhead, Improved particle filter
%     for nonlinear problems, IEE Proc.-Radar, Sonar Navig., Vol 146, No.
%     1, Februray 1999, pp 2-7
    
    % Fearnhead Appendix B (pp 898)
    N = length(weights);
    sampleidxs = logical(zeros(size(weights)));
    n_subset = sum(subsetidxs);
    
    if n_particles > n_subset
        error('srs:subsetTooSmall', 'n_particles must be no more than sum(subsetidxs)');
    end
    
    if length(subsetidxs) ~= length(weights)
        error('srs:wrongDims', 'length(subsetidxs) must equal length(weights)');
    end
    
    mean_target_weight = sum(weights(subsetidxs)) / n_particles;
    u = mean_target_weight * rand;    
    
    inds = find(subsetidxs);
    for ii = 1:n_subset
        i = inds(ii);
        u = u - weights(i);        
        if u < 0
            % Resample particle i means put it into our idxs
            sampleidxs(i) = true;
            u = u + mean_target_weight;
        end
    end   

end

