function [ blocks ] = block_sample_order( z, K, blocksz )
    N = length(z);
    Nblocks = ceil(length(z) / blocksz);
    blocks = zeros(Nblocks, blocksz + 1);
    capacity = zeros(Nblocks);
    
    idxs = randperm(N);
    
    
    for n = 1:N
        
    end
    

end

