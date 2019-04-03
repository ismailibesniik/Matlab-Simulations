%% MPC Building Control Project
%By Alex Karpilow and Besnik Ismaili



clc;
close all;


yalmip('clear')
clear all

%% Model data

load building.mat;
load battery.mat;

% Parameters of the Building Model
A = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C = ssM.C;
Ts = ssM.timestep;

% Parameters of the Storage Model
a = ssModel.A;
b = ssModel.Bu;   

% Installation Test
yalmip('version')
sprintf('The Project files are successfully installed')

%Plot disturbances:
figure;
plot(refDist(1,:));hold on; plot(refDist(2,:));plot(refDist(3,:));
legend('outside temp(C)','solar gains(kW)','internal gains (kW)')
xlabel('samples (20min)')

%% Controller Design (Setting-up MPC optimizer)
N=72;%horizon length
%Other parameters
T=576-N+1; %Simulation time
umax = 15; %kW
umin = 0;  %kW
ymax = 26; %C
ymin = 22; %C
nx   = size(A,1);
nu   = size(Bu,2);
ny   = size(C,1);
nd   = size(Bd,2);
u=sdpvar(nu,N-1,'full'); 
y=sdpvar(ny,N,'full');
d=sdpvar(nd,N,'full');
x=sdpvar(nx,N,'full');
% 
%% Section 1: tracking MPC %fill in here
yref = [24 24 24]';
R = 10;
obj=0;
con=[];
for i=1:N-1
    obj=obj+(y(:,i)-yref)'*R*(y(:,i)-yref);
    con=[con,...
        umin<=u(:,i)<=umax,...
        x(:,i+1)==A*x(:,i)+Bu*u(:,i)+Bd*d(:,i); %sys evolution
        y(:,i+1)==C*x(:,i+1)];
        con=[con,ymin<=y(:,i+1)<=ymax];
end
ops = sdpsettings('verbose',1);
controller1 = optimizer(con,obj,ops,[x(:,1);d(:)],u);

[xt1, yt1, ut1, t1] = simBuild(controller1, T, @shiftPred, N, 1);

%% Section 2: economic MPC and soft constraints
S = 10;  %penalty
c0 = 0.2; %fixed cost
obj = 0;
con = [];
eps = sdpvar(3,1,'full'); %single epsilon for all horizon steps
con=[eps>=0];
 
for i=1:N-1
    con=[con,...
        umin<=u(:,i)<=umax,...
         y(:,i+1)==C*x(:,i+1),...
         x(:,i+1)==A*x(:,i)+Bu*u(:,i)+Bd*d(:,i)]; %sys evolution
         con=[con,ymin-eps<=y(:,i+1)<=ymax+eps];
     obj = obj+sum(c0*u(:,i)) + eps'*S*eps; %quadratic penalty

end
ops = sdpsettings('verbose',1);
controller2 = optimizer(con,obj,ops,[x(:,1);d(:)],u);
[xt2, yt2, ut2, t2,cost_Econ] = simBuild(controller2, T, @shiftPred, N, 1);

%% Section 3: economic, soft constraints, and variable cost
S=10; %penalty
obj=0;
con=[];
eps=sdpvar(3,1,'full');
ct=sdpvar(N,1,'full');%variable cost

for i=1:N-1
    con=[con,...
        umin<=u(:,i)<=umax,...
        x(:,i+1)==A*x(:,i)+Bu*u(:,i)+Bd*d(:,i); %sys evolution
        y(:,i+1)==C*x(:,i+1)];
        con=[con, ymin-eps<=y(:,i)<=ymax+eps, eps>=0];
     obj = obj+sum(ct(i)*u(:,i))+eps'*S*eps;
end
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller3 = optimizer(con,obj,ops,[x(:,1);d(:);ct(:)],u);
[xt3, yt3, ut3, t3,cost_varCost] = simBuild(controller3, T, @shiftPred, N, 2);

%% Section 4 : Night setbacks
S=10; %penalty
obj=0;
con=[eps>=0];
eps=sdpvar(3,1,'full');
c_t=sdpvar(N,1,'full');%cost
sb=sdpvar(N,1,'full');%setbacks

for i=1:N-1
    con=[con,...
        umin<=u(:,i)<=umax,...
        x(:,i+1)==A*x(:,i)+Bu*u(:,i)+Bd*d(:,i); %sys evolution
        y(:,i+1)==C*x(:,i+1)];
        con=[con,ymin-eps-sb(i)<=y(:,i+1)<=ymax+eps+sb(i)];
     obj = obj+sum(c_t(i)*u(:,i))+eps'*S*eps;

end
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller4 = optimizer(con,obj,ops,[x(:,1);d(:);c_t(:);sb(:)],u);
[xt4, yt4, ut4, t4,cost_nightSB] = simBuild(controller4, T, @shiftPred, N, 3);

%% Section 5 : Battery coupled with the building
xb=sdpvar(N,1,'full');%state of charge
e=sdpvar(1,N-1,'full'); %power consumption from grid
v=sdpvar(1,N-1,'full'); %dis/charging power: e-sum(u)
c_t=sdpvar(N,1,'full');%cost
sb=sdpvar(N,1,'full');%night setbacks
%limits
xbmax=20;
xbmin=0;
vmax=20;
vmin=-20;
S=10; %penalty on constraint violation
obj=0;
con=[];

for i=1:N-1
    con=[con,...
        umin<=u(:,i)<=umax,...              %zone power input constraints
        xbmin<=xb(i+1)<=xbmax,...             %battery state constraints
        e(i)>=0,...                         %no injection into the grid
        vmin<=v(i)<=vmax,...                %dis/charging rate limits
        y(:,i+1)==C*x(:,i+1),...                %measured room temperatures
        v(i)==e(i)-sum(u(:,i)),...          %battery charging rate=pwr from grid -pwr consumption
        x(:,i+1)== A*x(:,i)+Bu*u(:,i)+Bd*d(:,i),... %sys evolution
        xb(i+1) == a*xb(i) + b*v(i)];       %battery state evolution
   con=[con, ymin-eps-sb(i)<=y(:,i+1)<=ymax+eps+sb(i), eps>=0];%zone temperature constraints
   obj = obj+c_t(i)*e(i)+eps'*S*eps;

end
ops = sdpsettings('verbose',1, 'solver', '+gurobi');
controller5 = optimizer(con,obj,ops,[x(:,1);xb(1);d(:);c_t(:);sb(:)],[u;v;e]);
[xt5, yt5, ut5, t5, et, xbt,cost_battery] = simBuildStorage(controller5, T, @shiftPred, N);
