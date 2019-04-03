classdef Building
    %Limitation caracteristics on the Heating System, Storage System,
    %Comfort constraints 
    properties
        umax = 15; %kW maximum power of the heating system per zone
        umin = 0;  %kW minimum power of the heating system per zone
        ymax = 26; %C upper comfort constraint
        ymin = 22; %C lower comfort constraint 
    end
end