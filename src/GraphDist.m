classdef GraphDist
% Compute geodesic distance. Separated from probability model code.
    
    properties
        misc_data
        
        % Hack; to softcode later
        x_max;
        y_max;
    end        
            
    methods
        function o = GraphDist(misc_data)
            o.misc_data = misc_data;
            o.x_max = size(misc_data.img, 1);
            o.y_max = size(misc_data.img, 2);
        end
        
        function ss = slic_segment(o, xy)
        % ss = get.slic_segment(xy) x(1) is x; x(2) is y coord
        
            s_x = round(xy(1));
            s_y = round(xy(2));
            
            % Clamp to the image boundary
            if s_x < 1 || s_x > o.x_max || s_y < 1 || s_y > o.y_max
                warning(['crossed image boundary: s_x = ' num2str(s_x) ...
                         ' s_y = ' num2str(s_y)]);
            end

            s_x = min(max(1, s_x), o.x_max);
            s_y = min(max(1, s_y), o.y_max);                

            ss = o.misc_data.segments(s_x,s_y);                        
        end
        
        function d = dist(o, X, mu)
        % d = dist(X, mu) return the precomputed geodesic distance from the
        % full-dimensional X and the cluster center mu
            sx = o.slic_segment(X(1:2));
            sm = o.slic_segment(mu(1:2));
            
            d = o.misc_data.slic_dist(sx, sm);            
        end
                        
        function d = slic_to_dist(o, x)                        
            if isempty(o.last_fitted_kss) || o.last_fitted_kss ~= o.k_slic_segment
                %warning('refit triggered');
                kidxs  = o.mdp.cluster_assigns == o.my_k; 
                kslics = o.mdp.X(kidxs,o.slic_dim);
                o.fit(kslics);
            end
            
            d = o.misc_data.slic_dist(o.last_fitted_kss,x);
        end
    end    
end

