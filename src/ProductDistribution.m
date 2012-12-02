classdef ProductDistribution < OnlineDistribution
% Applies independent probability models to distinct dimensions of X    
    
    properties
        idxranges;
        pdfs;
        
        N;
    end
    
    methods
        function o = ProductDistribution(varargin)
            % o = ProductDistribution(startidx, endidx, pdf, ...)
            % Construct a distribution that splits its input X into
            % components
            
            assert(mod(length(varargin), 3) == 0, 'varargin must be multiple of 3!');
            o.N = length(varargin) / 3;
            
            o.idxranges = cell(o.N, 1);
            o.pdfs      = cell(o.N, 1);            
            
            n = 1;
            k = 1;
            while k < length(varargin)
                o.idxranges{n} = varargin{k}:varargin{k+1};
                o.pdfs{n}      = varargin{k+2};                
                k = k + 3;
                n = n + 1;
            end
        end
        
        function fit(o, X)
            for i = 1:o.N                
                o.pdfs{i}.fit(X(:, o.idxranges{i}));
            end
        end
        
        function add_point(o, x)
            for i = 1:o.N
                o.pdfs{i}.fit(x(o.idxranges{i}));
            end
        end
        
        function remove_point(o, x)
            for i = 1:o.N
                o.pdfs{i}.remove_point(x(o.idxranges{i}));
            end
        end
        
        function p = pred_like(o, X)
            p = 1;
            for i = 1:o.N
                p = p * o.pdfs{i}.pred_like(X(:, o.idxranges{i}));
            end
        end
        
        function p = pred_like_scalar(o, x)
            p = 1;
            for i = 1:o.N
                p = p * o.pdfs{i}.pred_like_scalar(x(o.idxranges{i}));
            end
        end        
    end
    
end

