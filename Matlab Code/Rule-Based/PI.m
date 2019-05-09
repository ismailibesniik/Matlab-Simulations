function [U,integral] = PI(error,integral,umin,umax) 
    
Kp=[0.7;0.4;0.7];
Ki=[2;0.8;2];
Kd=0;

%PI controller
integral = integral + error;
U = Kp.*error + Ki.*integral; %+ Kd*derivative;

%saturation
for k=1:3
    if(U(k,:)<=umin)
        U(k,:) = umin;
    elseif(U(k,:)>=umax)
        U(k,:) = umax;
    end
end

end