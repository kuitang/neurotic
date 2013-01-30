classdef GraphDistDistribution < matlab.mixin.Copyable & OnlineDistribution
% Distribution of the geodesic (graph) distance of a given point to the
% cluster center 
% Applies independent probability models to distinct dimensions of X    
    
    properties
        underlying, mdp, misc_data;
        my_k;
        last_fitted_kss;
        
        slic_dim = 5;
        % Hack; to softcode later
        x_max = 200;
        y_max = 200;
    end
    
    properties (Dependent)
        % The slic segment OF THE CLUSTER CENTER.
        k_slic_segment;        
    end
    
    methods(Access = protected)                
        % Override copyElement method:
        function cpObj = copyElement(obj)            
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            cpObj.mdp = obj.mdp;
            cpObj.misc_data = obj.misc_data;
            cpObj.last_fitted_kss = obj.last_fitted_kss;
            
%            assert(~isempty(obj.my_k));
            cpObj.my_k = obj.my_k;
            
            % Make a deep copy of underlying
            cpObj.underlying = copy(obj.underlying);
        end
    end
    
    methods
        function kss = get.k_slic_segment(o)
            assert(o.my_k > 1);
            
            s_x = round(o.mdp.cluster_likes{o.my_k}.pdfs{1}.post_mean(1));
            s_y = round(o.mdp.cluster_likes{o.my_k}.pdfs{1}.post_mean(2));

            % Clamp to the image boundary
            if s_x < 1 || s_x > o.x_max || s_y < 1 || s_y > o.y_max
                warning(['crossed image boundary: s_x = ' num2str(s_x) ...
                         ' s_y = ' num2str(s_y)]);
            end

            s_x = min(max(1, s_x), o.x_max);
            s_y = min(max(1, s_y), o.y_max);                

            kss = o.misc_data.segments(s_x,s_y);            
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
        
        function o = GraphDistDistribution(rhs)
            if isfield(rhs, 'underlying') % copy constructor
%                assert(~isempty(rhs.my_k));
                o.underlying = copy(rhs.underlying);
                o.my_k = rhs.my_k;
                o.last_fitted_kss = rhs.last_fitted_kss;
            else % construction; rhs IS the underlying
                o.underlying = copy(rhs);
            end            
        end                
        
        function fit(o, x)
% fit(x) expects a vector x of slic segments. Will fit the *graph
% distances* of those slic segments to the slic segment of my_k.
            if isempty(x)
                o.underlying.fit([]);
            else            
                o.last_fitted_kss = o.k_slic_segment;
                d = o.misc_data.slic_dist(o.last_fitted_kss,x);
                o.underlying.fit(d);                                                
            end
        end
        
        function add_point(o, x)
            d = o.slic_to_dist(x);            
            o.underlying.add_point(d);
        end
        
        function remove_point(o, x)
            d = o.slic_to_dist(x);
            o.underlying.remove_point(d);            
        end
        
        function p = pred_like(o, x)
            d = o.slic_to_dist(x);
            % HACK: Left-censor
            % TODO: Read Kumar and Correct. (Maybe special-case the zero?)
            d = max(d, eps);
            p = o.underlying.pred_like(d);            
        end
        
        function p = pred_like_scalar(o, x)
            p = o.pred_like(x);            
        end
        
        function d = empirical_distribution(o)
            kidxs  = o.mdp.cluster_assigns == o.my_k; 
            kslics = o.mdp.X(kidxs,o.slic_dim);
            
            d = o.misc_data.slic_dist(o.k_slic_segment, kslics);
        end
        
    end
    
end

