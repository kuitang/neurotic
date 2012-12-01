classdef Uniform < OnlineDistribution
% Stupid dummy class...    
        
    properties
        one_over_N;
    end
    
    methods
        
        function o = Uniform(N)
            o.one_over_N = 1 / N;
        end
    
        function remove_point(o, n)
        end
        
        function add_point(o, n)
        end
        
        function fit(o, X)
        end
        
        function p = pred_like_scalar(o, x)
            p = o.one_over_N;
        end
        
        function p = pred_like(o, X)
            p = o.one_over_N;
        end
    end
    
end

