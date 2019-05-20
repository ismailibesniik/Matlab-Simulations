%% MPC Building Control
%Puprose: Master Thesis Project
%Author: Besnik Ismaili
%MPC with PV included, Normal battery, EWH, Reference tracking
%Without night setbacks!

clc;
close all;
yalmip('clear')
clear all
%% Model data
load building.mat;
load battery.mat;
load PV_power.mat;
load EWH_parameters.mat;
load PV_power;
%Increasing the power production to 150%
PV_power = power_PV;
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

%Limits of the HVAC system and Battery and EWH
umax    = 5; %kW -> Maximum power of the HP - COP = 3;
umin    = 0;  %kW -> Minimum power of the HP
ymax    = 23; %C  -> Maximum temperature in all the zones
ymin    = 19; %C  -> Minimum temperature in all the zones
SOCmax  = 20; %kWh -> Battery capacity
SOCmin  = 0; %kWh -> Minimum SOC
p_bmax  = 10.5/3; %kW -> Ramp up power of the battery
p_bmin  = -20; %kW -> Ramp down power of the battery
uemax = 10; %kW -> Maximum power of the EWH
uemin = 0; %kW -> Maximum power of the EWH
Tmin = 50; %Celsius Degrees
Tmax = 70; %Celsius Degrees
COP  = 3;
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

%Choose the algorithm to run for the the Energy Self-Consumption Maximization
%algorithm = 1; % -> MPC
%algorithm = 2; % -> Rule-Based
algorithm = 3; % -> Rule-Based with Load Scheduling

switch (algorithm)
    
    case 1
%% Setting-up MPC optimizer
N    = 72; %horizon length 72*20/60 = 24[hours]
T    = 576-N+1; %Simulation time
yref = [22 22 22]'; % Reference temperature in all the zones
R    = 10; % The weight on the reference temperature tracking
S    = 100000; %Penalty on constraint violation
reference = 0; % Choose 1 to track a reference temperature given by yref, or 0 to be between the boundaries ymin and ymax
[xt, yt, ut, t, pt, SOC, cost_grid, p_bt, cpt, uet, Tempt] = MPC_optimizer(A,Bu,C,Bd,N,T,yref,R,umax,umin,ymax,ymin,SOCmax,SOCmin,p_bmax,p_bmin,uemax,uemin,Tmax,Tmin,reference,COP);
    case 2
N    = 72; %horizon length 72*20/60 = 24[hours]
T    = 505; %Simulation time
% %Generating the 20 minutes time steps for the week under simulation
t1   = datetime('19-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
t2   = datetime('28-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
time = t1:minutes(20):t2;
[xt, yt, ut, t, pt, SOC,cost_grid, p_bt, cpt, uet, Tempt] = Rule_Based_optimizer(A,Bu,C,Bd,N,T,time,umax,umin,ymax,ymin,SOCmax,SOCmin,p_bmax,p_bmin,uemax,uemin,Tmax,Tmin,COP);
    case 3
N    = 72; %horizon length 72*20/60 = 24[hours]
T    = 576-N+1; %Simulation time
yref = [22 22 22]'; % Reference temperature in all the zones
R    = 10; % The weight on the reference temperature tracking
S    = 100000; %Penalty on constraint violation   
% %Generating the 20 minutes time steps for the week under simulation
t1   = datetime('19-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
t2   = datetime('28-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
time = t1:minutes(20):t2;
[xt, yt, ut, t, pt, SOC, p_bt, cpt, uet, Tempt] = Rule_Based_Scheduling(A,Bu,C,Bd,N,S,T,time,umax,umin,ymax,ymin,SOCmax,SOCmin,p_bmax,p_bmin,uemax,uemin,Tmax,Tmin,COP);
end
%% Calculate the energy self - consumption
index = find(pt>0); %Indexes of the consumed energy from the grid
Consumed_grid = sum(pt(index))/3; %Energy consumed from the grid
index1 = find(pt<0); %Indexes of the injected energy to the grid 
Injected_grid = (-1)*sum(pt(index1))/3; %Energy injected to the grid
Total_energy_consumed = sum(sum(ut))/3 + sum(uet); %Total energy consumed
Self_consumption = Total_energy_consumed - Consumed_grid; %Self-consumed energy
PV_energy = sum(PV_power(1:505))/3; %Energy produced by the PVs
cost_feed_in = Injected_grid*0.0543; %Feed-in profit
total_electricity_cost = cost_grid - cost_feed_in; %Total electricity cost

%Pie plot 
p = pie3([Self_consumption;Injected_grid],[0 1]);
labels = {'Self-Consumed','Feed-in'};
legend(labels,'Location','southoutside','Orientation','horizontal')
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Self-Consumed: ';'Feed-in: '}; 
combinedtxt = strcat(txt,percentValues); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
title('Self-consumed energy vs Feed-in');

%
figure;
p = pie3([Self_consumption;Consumed_grid],[0 1]);
labels = {'Self-Consumption','Grid-Consumption'};
legend(labels,'Location','southoutside','Orientation','horizontal')
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Self-Cons: ';'Grid-Cons: '}; 
combinedtxt = strcat(txt,percentValues); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
colormap([0 1 0;      %// red
          1 0 0]);       %// green
title('Total consumed energy');
%% Plot
pp = p_bt;
indexx = find(p_bt<0);
indexx1 = find(p_bt>0);
% p_bt(indexx) = p_bt(indexx)*(-1);
p_bt(indexx) = 0;
figure
%  plot(t,p_bt+uet+sum(ut,1),t,p_bt' + double(PV_power(1:505)))
%  t=t./3;
 plot(t,pp+uet+sum(ut,1),t,double(PV_power(1:505)))
% xlabel('Hours');
% ylabel('Storage State - PV Power');

% figure
% subplot(2,1,1)
% plot(t,xbt(1,:),t,power_PV(1:505))
% xlabel('Hours');
% ylabel('Storage State - PV Power');

% figure
% subplot(2,1,1)
% plot(t,xbt(1,:),t,power_PV(1:505))
% xlabel('Hours');
% ylabel('Storage State - PV Power');
% 
% % figure
% subplot(2,1,2)
% plot(t,pt(1,:),t,power_PV(1:505))
% hold on
% plot(t,10*cpt(1,:),'r',t,power_PV(1:505))
% legend('Electrical Power Purchased','High/Low Price Time','PV power')
% xlabel('Hours');
% ylabel('Power purchased (kW) - PV power');
% 
% 
% figure
% plot(t,p_bt(1,:),t,power_PV(1:505))
% hold on
% plot(t,10*cpt(1,:),'r')
% legend('Power to storage Purchased (v)','High/Low Price Time','PV power')
% xlabel('Hours');
% ylabel('Power input to the Storage (kW) - PV power'); 
% 
% figure
% plot(t,uet(1,:))
% hold on
% legend('Power to EWH')
% figure
% plot(t,Tempt(1,:))
% legend('Temperature evolution of EWH')