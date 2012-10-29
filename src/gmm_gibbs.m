function [ samples ] = gmm_gibbs(X, gmm, niters, burn_in, sample_freq, recovery_file )
    assert(nargin >= 5, 'not enough arguments!');
    
    [N, D] = size(X);
    
    samples = [];

    gmm.n = zeros(gmm.K, 1);
    for n = 1:N
        k = gmm.z(n);
        gmm.n(k) = gmm.n(k) + 1
    end
    
    if nargin >=6 && ~isempty(recovery_file)
        load(recovery_file);
        first = iter + 1;        
    else
        first = 1;
    end
    
    %figure
    %hold on
        
    for i = first:niters
        tic
        z_old = gmm.z;
                
        gmm = gmm_gibbs_iter(gmm, X);
        change = mean(abs(z_old ~= gmm.z));
        
        itertime = toc;
        rate     = N / itertime;        
        
        fprintf(1, 'Iter %d: time = %0.2f; rate = %.2f; change = %0.3f; loglike = %06.2f\n', ...
                i, itertime, rate, change, gmm.loglike);     

        if i > burn_in && mod(i, sample_freq) == 0
            samples = [samples gmm];
            iter = i;
            save(['iter' num2str(i) '.mat']);
        end        
        %bar(reshape(gmm.mean', 1, numel(gmm.mean)));
        %drawnow

    end
    
end


