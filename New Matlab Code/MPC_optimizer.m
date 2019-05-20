%% MPC Building Control
%Puprose: Master Thesis Project
%Author: Besnik Ismaili
%MPC with PV included, Normal battery, EWH, Without night setbacks!
%Function used in Energy_Self_Consumption_Max.m to perform the MPC
%controller
function [xt5, yt5, ut5, t5, pt5, xbt5, cost_grid5, p_bt5, cpt5, uet5, Tempt5] = MPC_optimizer(A,Bu,C,Bd,N,T,yref,R,umax,umin,ymax,ymin,SOCmax,SOCmin,p_bmax,p_bmin,uemax,uemin,Tmax,Tmin,reference,COP)
load building.mat;
load battery.mat;
load PV_power;
load data_meteo_swiss_2018_2019.mat
load EWH_parameters.mat
% Parameters of the Building Model
A  = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C  = ssM.C;
Ts = ssM.timestep; %20 minutes timestep

% Parameters of the Storage Model
a     = ssModel.A;
b     = ssModel.Bu;   
b_ch  = 0.97; %battery charging efficiency coefficient
b_dch = 1/1.1; %battery discharging efficiency coefficient

%EWH parameters
w_k = 0; % Representing the water usage. For now the hot water is not used. We only heat the Water Tank.
Tout = 60; %Temperature of the water going out of the Storage
Tin = 20; %Temperature of the water going in the Storage

nx   = size(A,1);
nu   = size(Bu,2);
ny   = size(C,1);
nd   = size(Bd,2);
u    = sdpvar(nu,N-1,'full'); 
ue   = sdpvar(1,N-1,'full'); 
y    = sdpvar(ny,N,'full');
d    = sdpvar(nd,N,'full');
x    = sdpvar(nx,N,'full');
PV   = sdpvar(1,N,'full');
% 
%% MPC - Battery and PV production coupled with the building
%Declaring the variables
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
%Limits
S       = 100000; %Penalty on constraint violation
M       = 10000;
obj     = 0;
con     = [];
plt = 1;%Choose to plot in the simBuildStorage_EWH.m
%
%% Constraints 
for i=1:N-1
    con=[con,...
               umin        <=    u(:,i)  <= umax,... %Zone power input constraints
               SOCmin      <=  SOC(i+1)  <= SOCmax,... %Battery state constraints
                 0         <=  p_plus(i) <= M*s(i),... % p_plus representing power consumed from the grid
           - M*(1-s(i))    <= p_minus(i) <= 0,... % p_min representing power fed-in to the grid
              p_bmin       <=    p_b(i)  <= p_bmax,... %Battery discharging/charging rate limits
                 0         <=   p_ch(i)  <= p_bmax*s_ch(i),... 
        p_bmin*(1-s_ch(i)) <=  p_dch(i)  <= 0,... 
                 0         <=    s(i)    <= 1,...
                 0         <=   s_ch(i)  <= 1,...  
               uemin       <=   ue(:,i)  <= uemax,... %EWH power input limits
               Tmin        <=  Temp(i+1) <= Tmax,... % EWH temperature limits
        p(i)      == p_plus(i) + p_minus(i),...   %power injected or consumed from the grid
        y(:,i+1)  == C*x(:,i+1),... %measured room temperatures
        p_b(i)    == PV(i) - sum(u(:,i)) - ue(:,i) + p(i) ,... %battery charging rate = pwr from PV - pwr consumption +- pwr from grid 
        p_b(i)    == p_ch(i) + p_dch(i),... %
        x(:,i+1)  == A*x(:,i) + Bu*u(:,i)*COP + Bd*d(:,i),... %sys evolution
        SOC(i+1)  == a*SOC(i) + b_ch*p_ch(i)/3 + b_dch*p_dch(i)/3,... %battery state evolution
        Temp(i+1) == CO*(CO1*Temp(i) + CO2*d(1,i) - CO3*w_k*(Tout - Tin) + CO4*ue(:,i))]; %EWH model
   
% con=[con,...
%                umin        <=    u(:,i)  <= umax,... %Zone power input constraints
%                SOCmin      <=  SOC(i+1)  <= SOCmax,... %Battery state constraints
%                  0         <=  p_plus(i) <= M*s(i),... % p_plus representing power consumed from the grid
%            - M*(1-s(i))    <= p_minus(i) <= 0,... % p_min representing power fed-in to the grid
%                  0         <=   p_ch(i)  <= p_bmax,... %Battery discharging/charging rate limits
%                  0         <=    s(i)    <= 1,...
%                uemin       <=   ue(:,i)  <= uemax,... %EWH power input limits
%                Tmin        <=  Temp(i+1) <= Tmax,... % EWH temperature limits
%         p(i)      == p_plus(i) + p_minus(i),...   %power injected or consumed from the grid
%         y(:,i+1)  == C*x(:,i+1),... %measured room temperatures
%         p_ch(i)    == PV(i) - sum(u(:,i)) - ue(:,i) + p(i) ,... %battery charging rate = pwr from PV - pwr consumption +- pwr from grid 
%         x(:,i+1)  == A*x(:,i) + Bu*u(:,i)*COP + Bd*d(:,i),... %sys evolution
%         SOC(i+1)  == a*SOC(i) + b_ch*p_ch(i)/3,... %battery state evolution
%         Temp(i+1) == CO*(CO1*Temp(i) + CO2*d(1,i) - CO3*w_k*(Tout - Tin) + CO4*ue(:,i));%EWH model
%         ]; 
   con=[con, ymin - eps <= y(:,i+1) <= ymax + eps , eps>=0];%zone temperature constraints
   obj = obj + price(i)*p_plus(i) + eps'*S*eps;
   if(reference == 1)
       obj = obj + (y(:,i+1) - yref)'*R*(y(:,i+1) - yref);
   end
end
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller5 = optimizer(con,obj,ops,[x(:,1);SOC(1,:);d(:);price(:);sb(:);PV(:);Temp(1,:)],[u;p_b;p;ue]);
[xt5, yt5, ut5, t5, pt5, xbt5, cost_grid5, p_bt5, cpt5, uet5, Tempt5] = simBuildStorage_EWH(controller5, T, @shiftPred, N,w_k,plt,COP);
end