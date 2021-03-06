%INPUTS:
%  controller - controller (optimizer instance) - the parameter order suggested for the controller is:
%           [x ; xb; d_pred(:) ; cp(:) ; sb(:)]            for economic MPC, where d_pred is the prediction of disturbance input d, cp is the electricity price vector, and sb is the night-setback offset, over the prediction horizon          
%  T - simulation time (in time-steps)
%  fhandle - function handle to the shiftPred function. The input to this
%               function is the time-step and the outputs of this function are the
%               predictions of disturbance d, electricity price cp, and the night-setback offset over the prediction horizon for the given time step. For
%               more details, check the documentation of this function.
%  N - Prediction Horizon of your MPC controller
%  option - 1 for simulation without variable cost and night-setbacks, 2 for variable cost, but no night-setbacks, and 3 for both variable cost and night setbacks

%OUTPUTS:
% xt - state as a function of time
% yt - output as a function of time
% ut - input as a function of time
% t - time (time-steps)




function [ xt, yt, ut, t ,cost_grid, cpt] = simBuild( controller, T, fhandle, N, option)

%T=simulation duration in time samples
load building.mat;
load battery.mat;
% Parameters of the Building Model
A = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C = ssM.C;
Ts =60;% ssM.timestep;

% Parameters of the Storage Model
a = ssModel.A;
b = ssModel.Bu; 
yref=[22;22;22];
umax = 15; %kW
umin = 0;  %kW

x = x0red;

nx = length(A);
nu = size(Bu,2);
nd = size(Bd,2);
ny = size(C,1);


xt = zeros(nx,T);   % T is the length of simulation in samples
yt = zeros(ny,T);
ut = zeros(nu,T);
t = zeros(1,T);
sbt = zeros(1,T);
cpt = zeros(1,T);

%cost
c=0.2;
cost_grid=0;
derivative = 0;
last_error = 0;
error= 0;
integral = [0;0;0];
%% Simulating the system and the controller
for i = 1:T 
[d_pred, cp, sb] = fhandle(i, N);
if option == 1          % No night-setbacks and no variable cost (example)
[U, id] = controller{[x; d_pred(:)]};                   
elseif option == 2      % Variable cost, but no night-setbacks
[U, id] = controller{[x; d_pred(:); cp(:)]};          
cpt(:,i) = cp(1,1);
elseif option == 3      % Variable cost and night-setbacks
[U, id] = controller{[x; d_pred(:); cp(:); sb(:)]};  
cpt(:,i) = cp(1,1);
sbt(:,i) = sb(1,1);
elseif option == 4
cpt(:,i) = cp(1,1);
y=C*x;
last_error = error;
error = yref - y;
% derivative = error - last_error;
% integral = integral + error;
% U = Kp.*error + Ki.*integral + Kd*derivative;
% for k=1:3
%     if(U(k,:)<=umin)
%         U(k,:) = umin;
%     elseif(U(k,:)>=umax)
%         U(k,:) = umax;
%     end
% end
[U,integral1] = PI(error,integral,umin,umax); 

integral = integral1;
% U = (C*x >= yref).*3;
% U = U + (C*x <= yref).*umax;
end
xt(:,i) = x;
ut(:,i) = U(1:nu,1);
yt(:,i) = C*x;

%calculate total costs
if option==1
    cost_grid=cost_grid+sum(c*ut(:,i));
elseif option==2 || option==3
    cost_grid=cost_grid+sum(cpt(:,i)*ut(:,i));
end

cost_grid = cost_grid + sum(cpt(:,i)*ut(:,i))/3;

t(1,i) = i;

disp(['Iteration ' int2str(i)])
% yalmiperror(id)

x = A*x + Bu*U(1:nu,1) + Bd*d_pred(:,1);

end

%% Generating the Plots

% Converting time scale from time-step to hours
t = t./3;


figure
% subplot(2,3,1)
plot(t, yt(1,:))
hold on
% if option == 3
%     hold on
%     plot(t, 26+sbt(1,:),'r')
%     plot(t, 22-sbt(1,:),'r')
%     legend('Zone-1 Temperature', 'Temperature Constraints')
% end
xlabel('Hours');
ylabel('Temperature (C)');


% figure
% subplot(2,3,2)
plot(t, yt(2,:), 'k')
% if option == 3
%     hold on
%     plot(t, 26+sbt(1,:),'r')
%     plot(t, 22-sbt(1,:),'r')
%     legend('Zone-2 Temperature', 'Temperature Constraints')
% end
xlabel('Hours');
ylabel('Temperature (C)');

% figure
% subplot(2,3,3)
plot(t, yt(3,:), 'c')
if option == 3
    hold on
    plot(t, 26+sbt(1,:),'r')
    plot(t, 22-sbt(1,:),'r')
    legend('Zone-1','Zone-2','Zone-3', 'Temperature Constraints')
else
    legend('Zone-1','Zone-2','Zone-3')
end
xlabel('Hours');
ylabel('Temperature (C)');


figure
% subplot(2,3,4)
plot(t,ut(1,:))
hold on
% if option == 2 || option == 3
%     hold on
%     plot(t,10*cpt(1,:),'r')
%     legend('Zone-1 Input', 'High/Low Price Time')
% end
xlabel('Hours');
ylabel('Power Input (kW)');


% figure
% subplot(2,3,5)
plot(t,ut(2,:),'k')
% if option == 2 || option == 3
%     hold on
%     plot(t,10*cpt(1,:),'r')
%     legend('Zone-2 Input', 'High/Low Price Time')
% end
xlabel('Hours');
ylabel('Power Input (kW)');


% figure
% subplot(2,3,6)
plot(t,ut(3,:),'c')
if option == 2 || option == 3
    hold on
    plot(t,10*cpt(1,:),'r')
    legend('Zone-1','Zone-2','Zone-3','High/Low Price Time')
else
    legend('Zone-1','Zone-2','Zone-3')
end
xlabel('Hours');
ylabel('Power Input (kW)');


end

