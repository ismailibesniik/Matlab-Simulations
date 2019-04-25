%% MPC Building Control
%Puprose: Master Thesis Project
%Author: Besnik Ismaili

clc;
close all;
yalmip('clear')
clear all

%% Model data

load building.mat;
load battery.mat;
load PV_power;
power_PV = 1.5*power_PV;
% Parameters of the Building Model
A  = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C  = ssM.C;
Ts = 60;%ssM.timestep; %20 minutes timestep

% Parameters of the Storage Model
a = ssModel.A;
b = ssModel.Bu;   

%Plot disturbances:
figure;
plot(refDist(1,:));hold on; plot(refDist(2,:));plot(refDist(3,:));hold on;plot(power_PV);
legend('outside temp(C)','solar gains(kW)','internal gains (kW)', 'PV power (kW)')
xlabel('samples (20min)')

figure;
plot(power_PV);
legend('PV power(kW)')
xlabel('samples (20min)')
% figure;
% plot(refDist(1,:));hold on; plot(refDist(2,:));plot(refDist(3,:));
% legend('outside temp(C)','solar gains(kW)','internal gains (kW)')
% xlabel('samples (20min)')
%% Controller Design (Setting-up MPC optimizer)
% 
N    = 72;%horizon length
T    = 576-N+1; %Simulation time
umax = 15; %kW
umin = 0;  %kW
ymax = 26; %C
ymin = 22; %C
yref = [22 22 22]';
R = 10;
nx   = size(A,1);
nu   = size(Bu,2);
ny   = size(C,1);
nd   = size(Bd,2);
u    = sdpvar(nu,N-1,'full'); 
y    = sdpvar(ny,N,'full');
d    = sdpvar(nd,N,'full');
x    = sdpvar(nx,N,'full');
PV   = sdpvar(1,N,'full');
% 
%% Convetional On/Off controller


[xt1, yt1, ut1, t1] = simBuild([], 505, @shiftPred, N, 4);

Total_energy_consumed = sum(sum(ut1))/3; %kwh
et = double(power_PV(1:505)') - sum(ut1);
index = find(et>0);
Injected_grid = sum(et(index))/3;
index1 = find(et<0);
Consumed_grid = (-1)*sum(et(index1))/3;
Self_consumption = Total_energy_consumed - Consumed_grid;
PV_energy = sum(power_PV(1:505))/3;


 