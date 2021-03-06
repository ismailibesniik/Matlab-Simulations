%% Rule-Based Algorithm for Self-Consumption Maximization and Electricity Bill Reduction
%Puprose: Master Thesis Project
%Author: Besnik Ismaili

clc;
close all;
yalmip('clear')
clear all
%% Model and PV data
load building.mat;
load battery.mat;
load PV_power;
load data_meteo_swiss_2018_2019.mat;
load EWH_parameters.mat;
meteo = str; % Create Column
%Converting the time on the structure of the data from meteo_swiss to datetime format
time_string = num2str(cell2mat(meteo.time));
time_string = [str2num(time_string(:,1:4)),str2num(time_string(:,5:6)), str2num(time_string(:,7:8)), str2num(time_string(:,9:10)), str2num(time_string(:,11:12))];
time_string(:,6) = 0;
meteo.datetime = datetime(time_string,'TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');

%Increasing the power production to 150%
PV_pred = power_PV';
PV_pred = double(PV_pred);

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

xb    = 0; %Initial battery state of charge
N     = 72;

nx = length(A);
nu = size(Bu,2);
nd = size(Bd,2);
ny = size(C,1);

%Generating the 20 minutes time steps for the week under simulation
t1   = datetime('19-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
t2   = datetime('28-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
time = t1:minutes(20):t2;
T    = 72;%505; %Number of samples

%Find indexes of global irradiation on the above mentioned dates
index1 = find(contains(str.time,'201902190000'));
index2 = find(contains(str.time,'201902280000'));
glob_irr_val = str2num(char(meteo.glob_irr)); %Global irradiation for the given date
glob_irr_coef = 17.13074617; %Global irradiation coefficient to PV prediciton
PV_pred_ann = glob_irr_coef*glob_irr_val; %Predicted PV production by ANN


xt    = zeros(nx,T);
yt    = zeros(ny,T);
ut    = zeros(nu,T);
t     = zeros(1,T);
SOC   = zeros(1,T);
Tempt = zeros(1,T); %EWH temperature
Tref  = zeros(1,T); %EWH temperature reference
uet   = zeros(1,T); %EWH input

et  = zeros(1,T);
xbt = zeros(1,T);
p_ch= zeros(1,T);

cpt  = zeros(1,T);
sbt  = zeros(1,T);
cost = 0;
xt(:,1) = x0red; %Initial state values
Tempt(1,1) = 50;
Tref(1,1) = 50;
ymax = [24;24;24];
ymin = [18;18;18];
umax = [15;15;15]; %kW -> Maximum power of the HP
umin = [0;0;0];  %kW -> Minimum power of the HP
Tmin = 50;
Tmax = 70;
uetmax = 9;
uetmin = 0;
SOC(1) = 0;
%% The main loop running the algorithm
for i = 1:T

    [~, cp, sb, ~] = shiftPred(i, N);
    d_pred = refDist;
    %Getting the date out of the time vector
    date = datestr(time(i));
    date = convertCharsToStrings(date);
    date = extractBefore(date,"-Feb-2019");
    date = str2num(date);
    
    %Creating the time intervals
    time22b = sprintf('%d-Feb-2019 22:00:00',date-1);
    time04  = sprintf('%d-Feb-2019 04:00:00',date);
    time06  = sprintf('%d-Feb-2019 06:00:00',date);
    time18  = sprintf('%d-Feb-2019 18:00:00',date);
    time22  = sprintf('%d-Feb-2019 22:00:00',date);
    time04n = sprintf('%d-Feb-2019 04:00:00',date+1);
    
    if(isbetween(time(i),time22b,time04)||isbetween(time(i),time22,time04n))
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
            yref = ymin; %Celsius degrees -> Reference to MIN
            Tref = Tmin;
%------ %Time between 04:00 - 06:00 ---------------------------------------
        case 2
            
            %Building temperature reference
            yref = ymax; %Celsius degrees -> Reference to MAX         
            Tref = Tmax; 
%------ %Time between 06:00 - 18:00 ---------------------------------------
        case 3
            
            %Check the PV energy production forecast for the next 20 min
            index = find(time(i) == meteo.datetime);
            Energy_PV_forecast = sum(PV_pred_ann(index:index+1,1))/2/3;
            
            %Devide the energy equally to the three zones  
            u_fore = [Energy_PV_forecast;Energy_PV_forecast;Energy_PV_forecast]./3;   
           
            %Checking if it will pass the maximum allowed room temperature 
            xref = C\ymax;
            x_check = A*xt(:,i) + Bu*u_fore + Bd*d_pred(:,i);
            y_check = C*x_check;
            
            ut(:,i) = u_fore;
            
            exceed_max = find(y_check>=ymax);
            ex_max_neg = find(y_check<=ymax);
            exceed_min = find(y_check<=ymin);
            ex_min_neg = find(y_check<=ymin);
            if(sum(sum(exceed_max)) ~= 0)
                
                  x1 = xt(:,i);
                  y = sdpvar(3,1,'full');
                  u = sdpvar(3,1,'full');
                  x = sdpvar(10,1,'full');
                  d = d_pred(:,i);
                  
                  obj = (y - ymax)'*1*(y - ymax);
                  con = [x == A*x1 + Bu*u + Bd*d,...
                         y == C*x,...
                  u(ex_max_neg,1) == u_fore(ex_max_neg,1)
                  umin <= u <= umax];
                  ops = sdpsettings('verbose',0);
                  optimize(con,obj,ops);
                  ut(:,i) = value(u);
            end
            
            if(sum(sum(exceed_min)) ~= 0)
                
                  x1 = xt(:,i);
                  y = sdpvar(3,1,'full');
                  u = sdpvar(3,1,'full');
                  x = sdpvar(10,1,'full');
                  d = d_pred(:,i);
                  
                  obj = (y - ymin)'*1*(y - ymin);
                  con = [x == A*x1 + Bu*u + Bd*d,...
                         y == C*x,...
                  u(ex_min_neg,1) == u_fore(ex_min_neg,1)
                  [0;0;0] <= u <= [15;15;15]];
                  ops = sdpsettings('verbose',0);
                  optimize(con,obj,ops);
                  ut(:,i) = value(u);
            end
            
            %Power left for the EWH
            if(sum(sum(ut(:,i))) < Energy_PV_forecast)    
                 uet(:,i) =  Energy_PV_forecast - sum(ut(:,i),1);
            end
            
            Tref(:,i+1) = 70; % Celsius degrees 
            ueref = (Tref(:,i+1)/CO - CO1*Tempt(:,i) - CO2*Troom + CO3)/CO4;
            
            ex_w_max = uet(:,i)>=ueref;
                
            if(ex_w_max ~= 0)
                uet(:,i) = ueref;  
            end
            
            Tref(:,i+1) = 50; % Celsius degrees 
            ueref = (Tref(:,i+1)/CO - CO1*Tempt(:,i) - CO2*Troom + CO3)/CO4;
            
            ex_w_min = uet(:,i)<=ueref;
                
            if(ex_w_min ~= 0)
                uet(:,i) = ueref;  
            end
            
%------ %Time between 18:00 - 22:00 ---------------------------------------
        case 4
            
            %Building temperature reference
            yref = ymin; %Celsius degrees
            Tref = Tmin;
    end
    
    if(mode == 1||mode == 2||mode==4)
        x1 = xt(:,i);
        y = sdpvar(3,1,'full');
        u = sdpvar(3,1,'full');
        x = sdpvar(10,1,'full');
        d = d_pred(:,i);
    %Calculating the input to track the reference -> HP
    obj = (y - yref)'*1*(y - yref);
    con = [x == A*x1 + Bu*u + Bd*d,...
           y == C*x,...
           [0;0;0] <= u <= [15;15;15]];
    ops = sdpsettings('verbose',0);
    optimize(con,obj,ops);
    ut(:,i) = value(u);
    
    %EWH model
    uet(:,i) = (Tref/CO - CO1*Tempt(:,i) - CO2*Troom + CO3)/CO4;
    
    p_ch(i)=0;
    end
    
    if ( mode == 3)
     %Limits on the input HP
    index11 = find(ut(:,i) < 0);
    index22 = find(ut(:,i) > 15);
    if (sum(sum(index11)))
       ut(index11,i) = 0;
    elseif(sum(sum(index22)))
       ut(index22,i) = 15;
    end
    end
    
    %Limits on the input EWH
    if (uet(:,i) < uetmin)
       uet(:,i) = uetmin;
    elseif(uet(:,i) > uetmax)
       uet(:,i) = uetmax;
    end
    
    if ( mode == 3)
    %Power left for the EV
    if(sum(sum(ut(:,i))) + uet(:,i) < Energy_PV_forecast)    
           p_ch(i) =  Energy_PV_forecast - sum(ut(:,i),1) - uet(:,i);
           
           if(p_ch(i) > 10.5)
               p_ch(i) = 10.5;
           end
    end
    end
    
    xt(:,i+1) = A*xt(:,i) + Bu*ut(:,i) + Bd*d_pred(:,i); 
    
    Tempt(:,i+1) = CO*(CO1*Tempt(:,i) + CO2*Troom - CO3 + CO4*uet(:,i));
    
    SOC(i+1)  = a*SOC(i) + b_ch*p_ch(i)/3;
    
    if(SOC(i+1)>100)
        SOC(i+1) = 100;
    end
    %Energy exchanged with the grid
    p(i) = PV_pred(:,i) - sum(ut(:,i),1);
    
    if(p(i)>0)
    feed_in(i) = p(i);%Feed-in energy
    else
    grid_cons(i) = (-1)*p(i);%Energy withdrawn from the grid
    end
    
    cpt(:,i) = cp(1,1);
    yt(:,i)  = C*xt(:,i);
    t(1,i)   = i;
    
    disp(['Iteration ' int2str(i)]);

end
t = t./3;
figure
% subplot(2,3,1)
figure;
plot(t, yt(:,:))
figure;
plot(t, Tempt(:,1:end-1))


