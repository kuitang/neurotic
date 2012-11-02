function [ h ] = plot_intensity_hists( mm, X )
% Plot intensity histogram for each class.
    h = mm.h_diagnostic;    
    
    for k = 1:mm.K        
        subplot(2,mm.K,k);
        idxs = mm.s_z == k;
        x = X(idxs,3);
        hist(x);        
        title([' N = ' num2str(sum(idxs)) ...
               ' m = ' num2str(mean(x), 2) ...
               ' sd = ' num2str(std(x, 1), 2)]);   
        
        % GMM only. TODO: Generalize the parametric plotter.
        if k > 1
            subplot(2,mm.K,mm.K + k);        
            ezplot(@(x) normpdf(x, mm.mean(k,3), sqrt(mm.cov(3,3,k))), linspace(0,1));
            title(['fitted m = ' num2str(mm.mean(k,3), 2) ...
                   ' sd = ' num2str(sqrt(mm.cov(3,3,k)), 2)]);
        end
        
        mm.cov
        
    end
        
end

