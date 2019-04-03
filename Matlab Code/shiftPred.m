%INPUT:
%  t - simulation time-step
%  N - Prediction Horizon of your MPC controller

%OUTPUTS:
% pred - prediction of the disturbance input, over the MPC prediction horizon
% cp - shifted price of electricity consumption over the MPC prediction horizon, at time step t
% sb - shifted comfort constraint off-sets over the MPC prediction horizon, at time step t


function [pred, cp, sb, PV_pred] = shiftPred(t, N)

load building.mat;
load PV_power;
%t is sample number (20min)

%% Disturbance Prediction
pred = refDist(:,t:t+N-1);
PV_pred = power_PV';
PV_pred = double(PV_pred(:,t:t+N-1));
%% Variable Price Prediction 
t_start=mod(t,72);%time step out of a day
if(t_start==0) %if exactly at start/end of day, set index=1
   t_start=1; 
end
t_end=t_start+N-1; %look 10 steps into future
cp1=[.2*ones(1,3*10),.04*ones(1,3*6), .2*ones(1,3*8)];%for 24 hours
cp1=[cp1 cp1]; %2 days
cp=cp1(t_start:t_end);

%% Night-Setback Prediction
%sdpvar instance representing the relaxation in the output over the
%predicition horizon
sb1=[4*ones(1,8*3),0*ones(1,10*3), 4*ones(1,6*3)];
sb1=[sb1 sb1]; %2 days
sb=sb1(t_start:t_end);



end

