classdef OnlineDistribution < matlab.mixin.Copyable
% OnlineDistribution supports add_point and remove_point functions    
    
    methods (Abstract)
        add_point(o, x)
        remove_point(o, x)
    end
    
end

