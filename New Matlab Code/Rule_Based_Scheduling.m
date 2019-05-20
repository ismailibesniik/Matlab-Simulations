%% Rule-Based Algorithm for Self-Consumption Maximization and Electricity Bill Reduction
%Puprose: Master Thesis Project
%Author: Besnik Ismaili
%Rule-Based Algorithm combined with Load Scheduling
function [xt, yt, ut, t, pt, xbt, p_bt, cpt, uet, Tempt] = Rule_Based_Scheduling(A,Bu,C,Bd,N,S,T,time,umax,umin,ymax,ymin,SOCmax,SOCmin,p_bmax,p_bmin,uemax,uemin,Tmax,Tmin,COP)
%% Model and PV data
load building.mat;
load battery.mat;
load PV_power;
load data_meteo_swiss_2018_2019.mat
load EWH_parameters.mat
meteo = str; % Create Column

% Parameters of the Building Model
A  = ssM.A;
Bu = ssM.Bu;
Bd = ssM.Bd;
C  = ssM.C;

% Parameters of the Storage Model
a     = ssModel.A;
b     = ssModel.Bu;   
b_ch  = 0.97; %battery charging efficiency coefficient
b_dch = 1/1.1; %battery discharging efficiency coefficient

%Converting the time on the structure of the data from meteo_swiss to datetime format
time_string = num2str(cell2mat(meteo.time));
time_string = [str2num(time_string(:,1:4)),str2num(time_string(:,5:6)), str2num(time_string(:,7:8)), str2num(time_string(:,9:10)), str2num(time_string(:,11:12))];
time_string(:,6) = 0;
meteo.datetime = datetime(time_string,'TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');

%Increasing the power production to 150%
PV_pred = power_PV';
PV_pred = double(PV_pred);

nx = length(A);
nu = size(Bu,2);
nd = size(Bd,2);
ny = size(C,1);
xt    = zeros(nx,T);
yt    = zeros(ny,T);
ut    = zeros(nu,T);
t     = zeros(1,T);
Tempt = zeros(1,T); %EWH temperature
Tref  = zeros(1,T); %EWH temperature reference
uet   = zeros(1,T); %EWH input

pt  = zeros(1,T);
p_ch = zeros(1,T);
xbt = zeros(1,T);
p_bt  = zeros(1,T);

cpt  = zeros(1,T);
sbt  = zeros(1,T);
cost = 0;
xt(:,1) = x0red; %Initial state values
SOC(1) = 0;  %Initial battery state of charge
Tempt(:,1) = 50;

%Find indexes of global irradiation on the above mentioned dates
glob_irr_val = str2num(char(meteo.glob_irr)); %Global irradiation for the given date
glob_irr_coef = 17.13074617; %Global irradiation coefficient to PV prediciton
PV_pred_ann = glob_irr_coef*glob_irr_val; %Predicted PV production by ANN


ymax = repmat(ymax,3,1); %Max temperature
ymin = repmat(ymin,3,1); %Min temperature
umax = repmat(umax,3,1); %kW -> Maximum power of the HP
umin = repmat(umin,3,1); %kW -> Minimum power of the HP
enter = 0;
k = 1;
M       = 10000;
%% The main loop running the algorithm
for i = 1:T

    %Getting the date out of the time vector
    date = datestr(time(i));
    date = convertCharsToStrings(date);
    date = extractBefore(date,"-Feb-2019");
    date = str2num(date);
    
    %Creating the time intervals
    time22b = sprintf('%d-Feb-2019 22:00:00',date-1);
    time044 = sprintf('%d-Feb-2019 03:59:00',date);
    time04  = sprintf('%d-Feb-2019 04:00:00',date);
    time06  = sprintf('%d-Feb-2019 06:00:00',date);
    time18  = sprintf('%d-Feb-2019 18:00:00',date);
    time22  = sprintf('%d-Feb-2019 22:00:00',date);
    time04n = sprintf('%d-Feb-2019 03:59:00',date+1);
    
    if(isbetween(time(i),time22b,time044)||isbetween(time(i),time22,time04n))
        mode = 1;
    elseif(isbetween(time(i),time04,time06))
        mode = 2;
    elseif(isbetween(time(i),time06,time18))
        mode = 3;
    elseif(isbetween(time(i),time18,time22))
        mode = 4;
    end
    
    switch mode
        
%------ %Time between 22:00 - 04:00 ---------------------------------------
        case 1
            
              %Building temperature reference
              if(enter==0)
              [U,d_pred,cp] = Load_Scheduler(xt(:,i),M,A,Bu,C,Bd,N,S,umax,umin,ymax,ymin,SOCmax,SOCmin,p_bmax,p_bmin,uemax,uemin,Tmax,Tmin,COP);
              enter = 1;
              end
              ut(:,i) = U(1:3,k);
%------ %Time between 04:00 - 06:00 ---------------------------------------
        case 2         
              %Building temperature reference    
              ut(:,i) = U(1:3,k);
%------ %Time between 06:00 - 18:00 ---------------------------------------
        case 3
              ut(:,i) = U(1:3,k);            
%------ %Time between 18:00 - 22:00 ---------------------------------------
        case 4  
            %Building temperature reference
            ut(:,i) = U(1:3,k);
            enter = 0;
    end
    
    xt(:,i+1) = A*xt(:,i) + Bu*ut(:,i)*COP + Bd*d_pred(:,k); 
     
    %Energy exchanged with the grid
    pt(i) = PV_pred(:,i) - sum(ut(:,i),1) - uet(:,i) - p_ch(i);
    
    if(pt(i)>0)
    feed_in(i) = pt(i);%Feed-in energy
    else
    grid_cons(i) = (-1)*pt(i);%Energy withdrawn from the grid
    end
    
    cpt(:,i) = cp(1,1);
    yt(:,i)  = C*xt(:,i);
    t(1,i)   = i;
    
    disp(['Iteration ' int2str(i)]);
    
    k = k + 1;
    if(mod(i,71) == 0)
    k = 1;%time step out of a day
    end
end
