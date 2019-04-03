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
% Parameters of the Building Model
A  = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C  = ssM.C;
Ts = 60;%ssM.timestep; %20 minutes timestep

% Parameters of the Storage Model
a = ssModel.A;
b = ssModel.Bu;   
b_ch = 0.97;
b_dch= 1.1;
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
yref = [24 24 24]';
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
%% MPC - Battery and PV production coupled with the building
% 
xb    = sdpvar(N,1,'full'); %state of charge
e_plus= sdpvar(1,N-1,'full'); %power consumption from grid
e_minus= sdpvar(1,N-1,'full'); %power given to the grid
e = sdpvar(1,N-1,'full'); %overall power consumed and given from/to the grid
v     = sdpvar(1,N-1,'full'); %dis/charging power: e-sum(u)
v_ch = sdpvar(1,N-1,'full'); %charging power
v_dch = sdpvar(1,N-1,'full'); %discharging power
c_t   = sdpvar(N,1,'full'); %cost
sb    = sdpvar(N,1,'full'); %night setbacks
eps   = sdpvar(3,1,'full'); %single epsilon for all horizon steps
bx    = sdpvar(N,1,'full'); %surplus variable
bch   = sdpvar(N,1,'full'); %surplus variable
%limits
xbmax = 20; %kWh
xbmin = 0; %kWh
vmax  = 20; %kW
vmin  = -20; %kW
S     = 100000; %penalty on constraint violation
obj   = 0;
con   = [];
M     = 10000;
for i=1:N-1
    con=[con,...
        umin<=u(:,i)<=umax,...              %zone power input constraints
        xbmin<=xb(i+1)<=xbmax,...           %battery state constraints
        0<=e_plus(i)<=M*bx(i),...           %no injection into the grid
        - M*(1-bx(i))<=e_minus(i)<=0,...    %surplus variable, injecting power from PV to the Grid
        e(i) == e_plus(i) + e_minus(i),...  %power injected or consumed from the grid
        vmin<=v(i)<=vmax,...                %dis/charging rate limits
        y(:,i+1)==C*x(:,i+1),...            %measured room temperatures
        v(i)==PV(i)-sum(u(:,i)) + e(i) ,... %battery charging rate=pwr from grid -pwr consumption
        v(i)==v_ch(i)+v_dch(i),...                         %to only charge the battery
        0<=v_ch(i)<=vmax*bch(i),... 
        vmin*(1-bch(i))<=v_dch(i)<=0,... 
        x(:,i+1)== A*x(:,i)+Bu*u(:,i)+Bd*d(:,i),... %sys evolution
        xb(i+1) == a*xb(i) + b_ch*v_ch(i)/3 + 1/b_dch*v_dch(i)/3,...%battery state evolution
        0<=bx(i)<=1,...
        0<=bch(i)<=1];     
   con=[con, ymin-eps-sb(i)<=y(:,i+1)<=ymax+eps+sb(i), eps>=0];%zone temperature constraints
   obj = obj+c_t(i)*e_plus(i)+eps'*S*eps;

end
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller5 = optimizer(con,obj,ops,[x(:,1);xb(1);d(:);c_t(:);sb(:);PV(:)],[u;v;e]);
[xt5, yt5, ut5, t5, et, xbt,cost_battery] = simBuildStorage(controller5, T, @shiftPred, N);


% [constraints, obj,u,y,d,x] = constraint_generator(A,Bu,C,Bd,N,yref,R,building);

%% Convetional On/Off controller


% [xt1, yt1, ut1, t1] = simBuild([], 505, @shiftPred, N, 4);
 