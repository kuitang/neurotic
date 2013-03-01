function [ cmat ] = threshCmat( blackdots, quantile )
    T = size(blackdots, 3);
    cmat = cell(T, 1);
    for t = 1:T
        [x y] = find(blackdots(:,:,t) > quantile);
        cmat{t} = [x y]
    end
end


