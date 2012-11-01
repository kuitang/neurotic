function [ h ] = plot_intensity_hists( mm, X )
% Plot intensity histogram for each class.
    h = mm.h_diagnostic;    
    
    for k = 1:mm.K        
        subplot(2,mm.K,k);
        idxs = mm.s_z == k;
        hist(X(idxs,3));
        title(['Class ' num2str(k) ' intensities']);   
        
        % GMM only. TODO: Generalize the parametric plotter.
        if k > 1
            subplot(2,mm.K,mm.K + k);        
            ezplot(@(x) normpdf(x, mm.mean(k,3), sqrt(mm.cov(3,3,k))), linspace(0,1));
            title(['Class ' num2str(k) ' fitted']);
        end
        
    end
        
end

