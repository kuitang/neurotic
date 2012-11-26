function [ cells ] = partition_rem( X, chunksz )
% Partition X by first dimension into into chunks of size chunksz, except for the remainder. 
% X must be at least a column vector.
    N = size(X, 1);
    assert(N > 0, 'cannot handle row vectors');
    n_chunks = floor(N / chunksz);
    r = rem(N, chunksz);
    
    parts = chunksz * ones(1, n_chunks);
    
    if r > 0
        n_chunks = n_chunks + 1;
        parts = [parts r];
    end
    
    cells = mat2cell(X, parts);

end

