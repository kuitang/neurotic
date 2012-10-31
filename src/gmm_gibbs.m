function [ samples ] = gmm_gibbs(X, img, gmm, niters, burn_in, sample_freq, recovery_file )
    assert(nargin >= 5, 'not enough arguments!');
    
    [N, D] = size(X);
    
    samples = [];

    gmm.n = zeros(gmm.K, 1);
    for n = 1:N
        k = gmm.s_z(n);
        gmm.n(k) = gmm.n(k) + 1;
    end
    
    if nargin >=7 && ~isempty(recovery_file)
        load(recovery_file);
        first = iter + 1;        
    else
        first = 1;
    end
        
    for i = first:niters
        tic
        z_old = gmm.s_z;
                
        gmm = gmm_gibbs_iter(gmm, X);
        change = mean(abs(z_old ~= gmm.s_z));
        
        itertime = toc;
        rate     = N / itertime;        
        
        fprintf(1, 'Iter %d: time = %0.2f; speed = %.2f; change = %0.3f; loglike = %06.2f\n', ...
                i, itertime, rate, change, gmm.loglike);     

        plot_point_overlay(img, gmm, X);
        drawnow;
        
        if i > burn_in && mod(i, sample_freq) == 0
            samples = [samples gmm];
            iter = i;
            if ~isempty(img)
                plot_point_overlay(img, gmm, X);
                drawnow;
            end
            %save(['iter' num2str(i) '.mat']);
        end        

    end
    
end


