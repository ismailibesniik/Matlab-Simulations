%Puprose: Master Thesis Project
%Author: Besnik Ismaili
%Load Scheduler function
function [U,d_pred,cp] = Load_Scheduler(xt,M,A,Bu,C,Bd,N,S,umax,umin,ymax,ymin,SOCmax,SOCmin,p_bmax,p_bmin,uemax,uemin,Tmax,Tmin,COP)
load building.mat;
load battery.mat;
load data_meteo_swiss_2018_2019.mat
load EWH_parameters.mat

% Parameters of the Storage Model
a     = ssModel.A;
b     = ssModel.Bu;   
b_ch  = 0.97; %battery charging efficiency coefficient
b_dch = 1/1.1; %battery discharging efficiency coefficient        
obj = 0;
con =[];
nx = length(A);
nu = size(Bu,2);
nd = size(Bd,2);
ny = size(C,1);
u    = sdpvar(nu,N-1,'full'); 
ue   = sdpvar(1,N-1,'full'); 
y    = sdpvar(ny,N,'full');
d    = sdpvar(nd,N,'full');
x    = sdpvar(nx,N,'full');
PV   = sdpvar(1,N,'full');
SOC     = sdpvar(N,1,'full'); %State of charge
p_plus  = sdpvar(1,N-1,'full'); %Power consumed from grid
p_minus = sdpvar(1,N-1,'full'); %Power fed in the grid
p       = sdpvar(1,N-1,'full'); %Power consumed from and fed in the grid
p_b     = sdpvar(1,N-1,'full'); %Battery dis/charging power: e-sum(u)
p_ch    = sdpvar(1,N-1,'full'); %Battery charging power
p_dch   = sdpvar(1,N-1,'full'); %Battery discharging power
price   = sdpvar(N,1,'full'); %Electricity price
sb      = sdpvar(N,1,'full'); %Night setbacks
eps     = sdpvar(3,1,'full'); %Single epsilon for all horizon steps
s       = sdpvar(N,1,'full'); %Surplus variable
s_ch    = sdpvar(N,1,'full'); %Surplus variable
Temp    = sdpvar(N,1,'full'); %EWH temperature
xb = 0; %initial battery state of charge
Temp1 = 50; %initial EWH temperature
w_k = 0;
%EWH parameters
Tin = 15; %-> [C degrees]
Tout = 60; %-> [C degrees]
for j=1:N-1
con=[con,...
               umin        <=    u(:,j)  <= umax,... %Zone power input constraints
               SOCmin      <=  SOC(j+1)  <= SOCmax,... %Battery state constraints
                 0         <=  p_plus(j) <= M*s(j),... % p_plus representing power consumed from the grid
           - M*(1-s(j))    <= p_minus(j) <= 0,... % p_min representing power fed-in to the grid
              p_bmin       <=    p_b(j)  <= p_bmax,... %Battery discharging/charging rate limits
                 0         <=   p_ch(j)  <= p_bmax*s_ch(j),... 
        p_bmin*(1-s_ch(j)) <=  p_dch(j)  <= 0,... 
                 0         <=    s(j)    <= 1,...
                 0         <=   s_ch(j)  <= 1,...  
               uemin       <=   ue(:,j)  <= uemax,... %EWH power input limits
               Tmin        <=  Temp(j+1) <= Tmax,... % EWH temperature limits
        p(j)      == p_plus(j) + p_minus(j),...   %power injected or consumed from the grid
        y(:,j+1)  == C*x(:,j+1),... %measured room temperatures
        p_b(j)    == PV(j) - sum(u(:,j)) - ue(:,j) + p(j) ,... %battery charging rate = pwr from PV - pwr consumption +- pwr from grid 
        p_b(j)    == p_ch(j) + p_dch(j),... %
        x(:,j+1)  == A*x(:,j) + Bu*u(:,j)*COP + Bd*d(:,j),... %sys evolution
        SOC(j+1)  == a*SOC(j) + b_ch*p_ch(j)/3 + b_dch*p_dch(j)/3,... %battery state evolution
        Temp(j+1) == CO*(CO1*Temp(j) + CO2*d(1,j) - CO3*w_k*(Tout - Tin) + CO4*ue(:,j))]; %EWH model
con=[con, ymin - eps <= y(:,j+1) <= ymax + eps, eps>=0];%zone temperature constraints
obj = obj + price(j)*p_plus(j) + eps'*S*eps;

end
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller = optimizer(con,obj,ops,[x(:,1);SOC(1,:);d(:);price(:);sb(:);PV(:);Temp(1,:)],[u;p_b;p;ue]);

[d_pred, cp, sb, PV_pred] = shiftPred(1, N);
[U, ~] = controller{[xt; xb; d_pred(:); cp(:); sb(:); PV_pred(:); Temp1]}; % this is the suggested form for the controller : you can change it provided buildSim.m is also accordingly changed

end