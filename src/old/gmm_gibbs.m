function [ samples ] = gmm_gibbs( varargin )   
    samples = run_mm(@gmm_gibbs_iter, varargin{:});
end

