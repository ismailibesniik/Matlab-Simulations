%A function checking if given input will pass the maximum allowed EWH temperature in one time step  
function [uet1] = max_temp_EWH(Tempt,uet,Tmax,Tmin,Tref,uetmin,uetmax,opt) 
            load EWH_parameters.mat;
            
            if (opt == 0)
            ueref = (Tmax/CO - CO1*Tempt - CO2*Troom + CO3)/CO4;
            uet1 = uet;
            ex_w_max = uet>=ueref;
                
            if(ex_w_max ~= 0)
                uet1 = ueref;  
            end
            
            ueref = (Tmin/CO - CO1*Tempt - CO2*Troom + CO3)/CO4;
            
            ex_w_min = uet<=ueref;
                
            if(ex_w_min ~= 0)
                uet1 = ueref;  
            end
            end
            
            if(opt == 1)
            uet1 = (Tref/CO - CO1*Tempt - CO2*Troom + CO3)/CO4;
            end
            
            %Check if exceeds the min or max allowed input of the system
            if (uet1 < uetmin)
            uet1 = uetmin;
            elseif(uet1 > uetmax)
            uet1 = uetmax;
            end
            
end