%INPUTS:
%  controller - controller (optimizer instance) - the parameter order suggested for the controller is:
%           [x ; xb; d_pred(:) ; cp(:) ; sb(:)]            for economic MPC, where d_pred is the prediction of disturbance input d, cp is the electricity price vector, and sb is the night-setback offset, over the prediction horizon          
%  T - simulation time (in time-steps)
%  fhandle - function handle to the shiftPred function. The input to this
%               function is the time-step and the outputs of this function are the
%               predictions of disturbance d, electricity price cp, and the night-setback offset over the prediction horizon for the given time step. For
%               more details, check the documentation of this function.
%  N - Prediction Horizon of your MPC controller
%

%OUTPUTS:
% xt - state as a function of time
% yt - output as a function of time
% ut - input as a function of time
% t - time (time-steps)
% et - electricity input to the battery model as a function of time
% xbt - State of the battery storage as a function of time


function [ xt, yt, ut, t, et, xbt, cost_battery, vt, cpt] = simBuildStorage(controller, T, fhandle, N)
load building.mat;
load battery.mat;
load PV_power;
% Parameters of the Building Model
A = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C = ssM.C;
Ts = ssM.timestep;

% Parameters of the Storage Model
a = ssModel.A;
b = ssModel.Bu; 
b_ch = 0.97;
b_dch= 1.1;

x = x0red;
xb = 0;%initial battery state of charge

nx = length(A);
nu = size(Bu,2);
nd = size(Bd,2);
ny = size(C,1);


xt = zeros(nx,T);   % T is the length of simulation in samples
yt = zeros(ny,T);

ut = zeros(nu,T);
t = zeros(1,T);


et = zeros(1,T);
xbt = zeros(1,T);
vt = zeros(1,T);

cpt = zeros(1,T);
sbt = zeros(1,T);
cost_battery=0;

for i = 1:T
[d_pred, cp, sb, PV_pred] = fhandle(i, N);
[U, id] = controller{[x; xb; d_pred(:); cp(:); sb(:); PV_pred(:)]}; % this is the suggested form for the controller : you can change it provided buildSim.m is also accordingly changed

xt(:,i) = x;
ut(:,i) = U(1:nu,1);
et(:,i) = U(end,1);
vt(:,i) = U(end-1,1);
xbt(:,i) = xb;
cpt(:,i) = cp(1,1);
sbt(:,i) = sb(1,1);

yt(:,i) = C*x;
t(1,i) = i;
cost_battery=cost_battery+cpt(:,i)*et(:,i);

disp(['Iteration ' int2str(i)]);
yalmiperror(id);
x = A*x + Bu*ut(:,i) + Bd*d_pred(:,1);
if(vt(:,i)>0)
xb = a*xb + b_ch*vt(:,i)/3;
else
xb = a*xb + 1/b_dch*vt(:,i)/3;
end

end

%% Generating the Plots

% Converting time scale from time-step to hours
 t = t./3;

figure
% subplot(2,3,1)
plot(t, yt(1,:))
hold on
% plot(t, 26+sbt(1,:),'r')
% plot(t, 22-sbt(1,:),'r')
% legend('Zone-1 Temperature', 'Temperature Constraints')
% xlabel('Hours');
% ylabel('Temperature - Zone1 (C)');


% figure
% subplot(2,3,2)
plot(t, yt(2,:),'k')
% hold on
% plot(t, 26+sbt(1,:),'r')
% plot(t, 22-sbt(1,:),'r')
% legend('Zone-2 Temperature', 'Temperature Constraints')
% xlabel('Hours');
% ylabel('Temperature(C)');

% figure
% subplot(2,3,3)
plot(t, yt(3,:),'c')
hold on
plot(t, 26+sbt(1,:),'r')
plot(t, 22-sbt(1,:),'r')
%%
% plot(t,eps_save(1,:));
% plot(t,eps_save(2,:));
% plot(t,eps_save(3,:));
% legend('Zone-3 Temperature', 'Temperature Constraints','viol zone1','viol zone2','viol zone3')
legend('Zone-1','Zone-2','Zone-3', 'Temperature Constraints', 'viol zone1','viol zone2','viol zone3')
xlabel('Hours');
ylabel('Temperature(C)');



figure
% subplot(2,3,4)
plot(t,ut(1,:))
hold on
% plot(t,10*cpt(1,:),'r')
% legend('Zone-1 Input', 'High/Low Price Time')
xlabel('Hours');
ylabel('Power Input - Zone1 (kW)');


% figure
% subplot(2,3,5)
plot(t,ut(2,:),'b')
% hold on
% plot(t,10*cpt(1,:),'r')
% legend('Zone-2 Input', 'High/Low Price Time')
xlabel('Hours');
ylabel('Power Input - Zone2 (kW)');


% figure
% subplot(2,3,6)
plot(t,ut(3,:),'c')
hold on
plot(t,10*cpt(1,:),'r')
% legend('Zone-3 Input', 'High/Low Price Time')
legend('Zone-1','Zone-2','Zone-3','High/Low Price Time')
xlabel('Hours');
ylabel('Power Input - Zone3 (kW)');



figure
subplot(2,1,1)
plot(t,xbt(1,:))
xlabel('Hours');
ylabel('Storage State');

% figure
subplot(2,1,2)
plot(t,et(1,:))
hold on
plot(t,10*cpt(1,:),'r')
legend('Electrical Power Purchased','High/Low Price Time')
xlabel('Hours');
ylabel('Power purchased (kW)');


figure
plot(t,vt(1,:))
hold on
plot(t,10*cpt(1,:),'r')
legend('Power to storage Purchased (v)','High/Low Price Time')
xlabel('Hours');
ylabel('Power input to the Storage (kW)');


end

