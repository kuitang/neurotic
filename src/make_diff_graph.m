function [ segments, slic_inds, edge_G ] = make_diff_graph( img, slic_npixels, slic_reg )
% make_diff_graph - returns sparse adjacency graph with pixel-difference edges
%
% img          - [X Y] image matrix
% slic_npixels - region size passed to vl_slic
% slic_reg     - regularizer passed to vl_slic
%
% segments - [X Y] matrix mapping pixels to superpixel ids
% edge_G   - sparse [Nslic Nslic] lower triangular graph. Two superpixels
%            are adjacent if there exist two adjacent underlying pixels.

    [R, C] = size(img);
    
    % vl_slic returns 0-indexed; we want 1-indexed.
    segments = vl_slic(single(img), slic_npixels, slic_reg) + 1;
    
%     figure(15);
%     imagesc(segments);
    
    mean_img = showslic(segments, img);
%     figure(16);
%     imshow(mean_img);
    
    assert(min(segments(:)) > 0);
    Nslic = max(segments(:));            
    
    % Allocate some extra to be on the safe side; doesn't really matter
    ivec = zeros(Nslic^2,1); jvec = zeros(Nslic^2,1); svec = zeros(Nslic^2,1);
    
    % Neighborhoods (ignore the edges)
    sz = size(img);
    
    nedge = 1;
    for n = 1:Nslic
        inds = uint32(find(segments(:) == n));
        slic_inds{n} = inds;
        
        neighbor_w = zeros(Nslic, 1);        
        
        for ind = inds'
            [r, c] = ind2sub(sz, ind);           
                
            % Explore the neighborhood
            for dr = [-1 1]
                for dc = [-1 1]
                    if (r + dr > 1 && r + dr < R) && (c + dc > 1 && c + dc < C)
                        other_n = segments(r + dr, c + dc);                    

                        % CHECK THIS LINE... earlier you asserted other_n >
                        % 0; should really check if its a different
                        % segment.
                        %
                        % Should really take difference between superpixel
                        % means.
                        if other_n ~= n
                            neighbor_w(other_n) = abs(mean_img(r + dr, c + dc) - mean_img(r, c));
                                                        
%                             if abs(mean_img(r + dr, c + dc) - mean_img(r, c)) > 0.05
%                                 neighbor_w(other_n) = 1;
%                             else
%                                 neighbor_w(other_n) = 0.01;
%                             end                                                        
                        end
                    end
                end
            end
        end

        % Now add edges                
        neighbor_idxs = find(neighbor_w);

        next_nedge = nedge + length(neighbor_idxs);
        ivec(nedge:(next_nedge-1)) = n;
        jvec(nedge:(next_nedge-1)) = neighbor_idxs;
        svec(nedge:(next_nedge-1)) = neighbor_w(neighbor_idxs);

        % Prepare for next iteration
        nedge = next_nedge;        
    end

    % Squeeze nonzeros
    nzidxs = ivec > 0;
    ivec = ivec(nzidxs);
    jvec = jvec(nzidxs);
    svec = svec(nzidxs);
    edge_G = sparse(ivec, jvec, svec, double(Nslic), double(Nslic)); 

end

