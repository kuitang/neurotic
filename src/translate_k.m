function [ k_idx, k_invidx ] = translate_k( gmm )
    nzidxs   = gmm.n > 0;
    assert(sum(nzidxs) == gmm.K, 'something is horribly wrong');
    k_idx    = find(nzidxs);
    k_invidx = zeros(size(gmm.n));    
    % relabel each true class label to a compact class label
    k_invidx(nzidxs) = 1:gmm.K;
end

