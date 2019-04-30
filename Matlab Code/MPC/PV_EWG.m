%% MPC Building Control
%Puprose: Master Thesis Project
%Author: Besnik Ismaili

% clc;
% close all;
% yalmip('clear')
% clear all
% %% Model data
% load building.mat;
% load battery.mat;
% load PV_power;
% load results
% 
% Npvs = 15; % - the number of PV modules in series
% Npvp = 6;% - the number of PV modules in paralel
% NOCT = 45.5; %-> [1/C degrees] - the nominal operating cell temperatur
% gamma = 0.00043;%-> [1/C degrees] - the power temperature coefficient at the maximum power point (MPP)
% Pv_stc = 165; %-> [W] - the power output of a PV-module under standard test conditions(STC)
% I_stc = 1000; %-> [W/m^2] - the irradiation at STC (Standard Testing Conditions)
% Tj_stc = 25; %-> [C degrees] - the cell temperature at STC
% I_noct = 800; %-> [W/m^2] - the irradiation at nominal operating temperature
% T_amb_NOCT = 20; %-> [C degrees] 
% %Pv(t) - total power output at time t of all the PV-modules.
% %I - the current global irradiation on the sloped surface of the modules. - the power temperature coefficient at the maximum power point (MPP)
% %Tamb(t) - the outside temperature.
% date = results.data(1).timestamp_local;
% results.data(1).timestamp_local
% I = zeros(49,1);
% T_amb = zeros(49,1);
% Tj = zeros(49,1);
% Pv = zeros(49,1);
% 
% for t = 1:1:49
% I(t) = results.data(t).solar_rad;
% T_amb(t)  = results.data(t).temp;
% %PV modeling 
% Tj(t) = T_amb(t) + I(t)/I_noct*(NOCT - T_amb_NOCT);
% Pv(t) = Pv_stc*I(t)/I_stc*(1-gamma*(Tj(t) - Tj_stc))*Npvs*Npvp;
% end
% 
% t1 = datetime('15-Apr-2019 10:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
% t2 = datetime('17-Apr-2019 10:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
% time = t1:hours(1):t2;
% 
% figure;
% plot(time, Pv);
% hold on
% plot(time,I);
% hold on
% plot(time,T_amb)
%% MPC Building Control
%Puprose: Master Thesis Project
%Author: Besnik Ismaili
%MPC with an Electrical Water Heater in the system.

clc;
close all;
yalmip('clear')
clear all
%% Model data
load building.mat;
load battery.mat;
load PV_power;

%Increasing the power production to 150%
PV_power = 1.5*power_PV; 

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

%Time of the year observed
t1 = datetime('20-Feb-2018','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
t2 = datetime('28-Feb-2018','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
t = (t1:minutes(20):t2)';
t = t(1:576); %preparing the time vector for plotting

%Plot disturbances:
figure;
plot(t,refDist(1,:));hold on; plot(t,refDist(2,:)/1000);plot(t,refDist(3,:));hold on;plot(t,power_PV);
legend('outside temp(C)','solar gains(kW)','internal gains (kW)', 'PV power (kW)')
xlabel('samples (20min)')
%% Controller Design (Setting-up MPC optimizer) 
N    = 72; %horizon length 72*20/60 = 24[hours]
T    = 576-N+1; %Simulation time
yref = [22 22 22]'; % -> Reference temperature in all the zones
R    = 10;
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
% w_k     = rand(N,1);
% w_k     = repmat(w_k,8,1);
w_k = 0; % Representing the water usage. For now the hot water is not used. We only heat the Water Tank.

%Limits
umax    = 15; %kW -> Maximum power of the HP
umin    = 0;  %kW -> Minimum power of the HP
ymax    = 24; %C  -> Maximum temperature in all the zones
ymin    = 20; %C  -> Minimum temperature in all the zones
SOCmax  = 20; %kWh -> Battery capacity
SOCmin  = 0; %kWh -> Minimum SOC
p_bmax  = 20; %kW -> Ramp up power of the battery
p_bmin  = -20; %kW -> Ramp down power of the battery
S       = 100000; %Penalty on constraint violation
obj     = 0;
con     = [];
M       = 10000;
uemax = 4.5; %kW -> Maximum power of the HP
uemin = 0; %kW -> Maximum power of the HP
Tmin = 60; %Celsius Degrees
Tmax = 80; %Celsius Degrees

%EWH parameters
a = 128.38; %-> [J/min C degrees]
c_w = 4.1813; %-> [J/g C degrees]
m_w = 196.82; %-> [kg]
C1 = 8.22*10^5; %-> [J/C degrees]
Tin = 10; %-> [C degrees]
Tout = 60; %-> [C degrees]
Troom = 22; %-> [C degrees]
Pmax = 4.5; %-> [kW]
%Constraints 
for i=1:N-1
    con=[con,...
               umin        <=  u(:,i)  <= umax,... %Zone power input constraints
               SOCmin      <=  SOC(i+1)  <= SOCmax,... %Battery state constraints
                 0         <=  p_plus(i) <= M*s(i),... % p_plus representing power consumed from the grid
           - M*(1-s(i))    <= p_minus(i) <= 0,... % p_min representing power fed-in to the grid
              p_bmin       <=    p_b(i)  <= p_bmax,... %Battery discharging/charging rate limits
                 0         <=   p_ch(i)  <= p_bmax*s_ch(i),... 
        p_bmin*(1-s_ch(i)) <=  p_dch(i)  <= 0,... 
                 0         <=      s(i)  <= 1,...
                 0         <=   s_ch(i)  <= 1,...  
               uemin       <=   ue(:,i)  <= uemax,...
               Tmin        <=  Temp(i+1) <= Tmax,...
        p(i)      == p_plus(i) + p_minus(i),...   %power injected or consumed from the grid
        y(:,i+1)  == C*x(:,i+1),... %measured room temperatures
        p_b(i)    == PV(i) - sum(u(:,i)) - ue(:,i) + p(i) ,... %battery charging rate = pwr from PV - pwr consumption +- pwr from grid 
        p_b(i)    == p_ch(i) + p_dch(i),... %to only charge the battery
        x(:,i+1)  == A*x(:,i) + Bu*u(:,i) + Bd*d(:,i),... %sys evolution
        SOC(i+1)  == a*SOC(i) + b_ch*p_ch(i)/3 + b_dch*p_dch(i)/3,... %battery state evolution
        Temp(i+1) == (1+a/(2*C1))^(-1)*((1-a/(2*C1))*Temp(i) + a/C1*Troom - 1/m_w*w_k*(Tout - Tin) + 20/C1*ue(i)*1000)%EWH model
        ];   
   con=[con, ymin - eps <= y(:,i+1) <= ymax + eps , eps>=0];%zone temperature constraints
   obj = obj + price(i)*p_plus(i) + eps'*S*eps + (y(:,i+1) - yref)'*R*(y(:,i+1) - yref);

end
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller5 = optimizer(con,obj,ops,[x(:,1);SOC(1,1);d(:);price(:);sb(:);PV(:);Temp(1,1)],[u;p_b;p;ue]);
[xt5, yt5, ut5, t5, et, xbt, cost_grid, vt5, cpt5, ue5, Temp5] = simBuildStorage_EWH(controller5, T, @shiftPred, N,w_k);
%% Calculate the energy self - consumption

index = find(et>0); %Indexes of the consumed energy from the grid
Consumed_grid = sum(et(index))/3; %Energy consumed from the grid
index1 = find(et<0); %Indexes of the injected energy to the grid 
Injected_grid = (-1)*sum(et(index1))/3; %Energy injected to the grid
Total_energy_consumed = sum(sum(ut5))/3; %Total energy consumed
Self_consumption = Total_energy_consumed - Consumed_grid; %Self-consumed energy
PV_energy = sum(PV_power(1:505))/3; %Energy produced by the PVs
cost_feed_in = Injected_grid*0.0543; %Feed-in profit
total_electricity_cost = cost_grid - cost_feed_in; %Total electricity cost
%% Plot
figure
subplot(2,1,1)
plot(t5,xbt(1,:),t5,power_PV(1:505))
xlabel('Hours');
ylabel('Storage State - PV Power');

% figure
subplot(2,1,2)
plot(t5,et(1,:),t5,power_PV(1:505))
hold on
plot(t5,10*cpt5(1,:),'r',t5,power_PV(1:505))
legend('Electrical Power Purchased','High/Low Price Time','PV power')
xlabel('Hours');
ylabel('Power purchased (kW) - PV power');


figure
plot(t5,vt5(1,:),t5,power_PV(1:505))
hold on
plot(t5,10*cpt5(1,:),'r')
legend('Power to storage Purchased (v)','High/Low Price Time','PV power')
xlabel('Hours');
ylabel('Power input to the Storage (kW) - PV power'); 

figure
plot(t5,ue5(1,:))
hold on

figure
plot(t5,Temp5(:,1))
