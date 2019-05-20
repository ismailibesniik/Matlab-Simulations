%Function to check if the given input in one time step will pass the max
%allowed temperature of the system
%The parameter optimize is optimize = 1 if we only want to find the optimal u to
%reach a reference temperature or optimize = 0 if we want to check if for a
%given input we exceed the max or min allowed temperature
function [u1] = max_temp_building(A,Bu,Bd,C,xt,ut,dt,ymax,ymin,yref,umax,umin,nx,ny,nu,opt,COP)
            

            y = sdpvar(ny,1,'full');
            u = sdpvar(nu,1,'full');
            x = sdpvar(nx,1,'full');
         if(opt == 0)
            u1 = ut;
            x_check = A*xt + Bu*ut*COP + Bd*dt;
            y_check = C*x_check;
            exceed_max = find(y_check>=ymax);
            ex_max_neg = find(y_check<=ymax);
            exceed_min = find(y_check<=ymin);
            ex_min_neg = find(y_check>=ymin);
            
            if(sum(sum(exceed_max)) ~= 0)
 
                  obj = (y - ymax)'*1*(y - ymax);
                  con = [x == A*xt + Bu*u*COP + Bd*dt,...
                         y == C*x,...
                         u(ex_max_neg,1) == ut(ex_max_neg,1),...
                         umin <= u <= umax];
                  ops = sdpsettings('verbose',0);
                  optimize(con,obj,ops);
                  u1 = value(u);
            end
            
            if(sum(sum(exceed_min)) ~= 0)
                  
                  obj = (y - ymin)'*1*(y - ymin);
                  con = [x == A*xt + Bu*u*COP + Bd*dt,...
                         y == C*x,...
                  u(ex_min_neg,1) == u1(ex_min_neg,1)
                  umin <= u <= umax];
                  ops = sdpsettings('verbose',0);
                  optimize(con,obj,ops);
                  u1 = value(u);
            end
         end
         
         if (opt == 1) 
             obj = (y - yref)'*1*(y - yref);
             con = [x == A*xt + Bu*u*COP + Bd*dt,...
                    y == C*x,...
                    umin <= u <= umax];
             ops = sdpsettings('verbose',0);
             optimize(con,obj,ops);
             u1 = value(u);
         end
         
         %Check if exceeds the min or max allowed input of the system
         index11 = find(u1 < 0);
         index22 = find(u1 > 5);
         if (sum(sum(index11)))
         u1(index11) = 0;
         elseif(sum(sum(index22)))
         u1(index22) = 5;
         end
end