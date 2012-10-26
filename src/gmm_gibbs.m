function [ samples ] = gmm_gibbs(X, gmm, niters, burn_in, sample_freq, recovery_file )
    assert(nargin >= 5, 'not enough arguments!');
    
    [N, D] = size(X);
    
    samples = [];

    if nargin >=6 && ~isempty(recovery_file)
        load(recovery_file);
        first = iter + 1;        
    else
        first = 1;
    end
    
    tic
    for i = first:niters
        tic
        z_old = gmm.z;
        gmm = gmm_gibbs_iter(gmm, X);
        change = mean(abs(z_old ~= gmm.z));
        disp(['Iteration ' num2str(i) ' took ' num2str(toc) ' seconds; ' num2str(change) ' changed classes']);        

        if i > burn_in && mod(i, sample_freq) == 0
            samples = [samples gmm];
            iter = i;
            save(['iter' num2str(i) '.mat']);
        end
        
        gmm.n
        gmm.mu        

    end
  
end


