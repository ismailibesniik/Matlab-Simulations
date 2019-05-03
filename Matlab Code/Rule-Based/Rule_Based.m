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

%Increasing the power production to 150%
PV_power = 1.5*power_PV; 

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

x     = x0red; %Initial state values
xb    = 0; %Initial battery state of charge
Temp1 = 70; %Initial EWH temperature
N     = 72;

nx = length(A);
nu = size(Bu,2);
nd = size(Bd,2);
ny = size(C,1);

%Generating the 20 minutes time steps for the week under simulation
t1   = datetime('19-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
t2   = datetime('28-Feb-2019 00:00:00','TimeZone','Europe/Zurich','Format','dd-MMM-yyyy HH:mm:ss');
time = t1:minutes(20):t2;
T    = size(time,2); %Number of samples

xref  = zeros(nx,T); % T is the length of simulation in samples
xt    = zeros(nx,T);
yt    = zeros(ny,T);
ut    = zeros(nu,T);
t     = zeros(1,T);
Tempt = zeros(1,T); %EWH temperature
Tref  = zeros(1,T); %EWH temperature
uet   = zeros(1,T); %EWH input

et  = zeros(1,T);
xbt = zeros(1,T);
vt  = zeros(1,T);

cpt  = zeros(1,T);
sbt  = zeros(1,T);
cost = 0;

%EWH parameters
a1    = 128.38; %-> [J/min C degrees]
c_w   = 4.1813; %-> [J/g C degrees]
m_w   = 196.82; %-> [kg]
C1    = 8.22*10^5; %-> [J/C degrees]
Tin   = 10; %-> [C degrees]
Tout  = 60; %-> [C degrees]
Pmax  = 4.5; %-> [kW]
Troom = 22; %-> [C degrees]
w_k = 0; %Water usage variable
CO = ((1+20*30*a1/(60*2*C1))^(-1));
CO1= (1-20*30*a1/(60*2*C1));
CO2= 20*30*a1/(60*C1);
CO3= 20*30*1/m_w*w_k*(Tout - Tin);
CO4= 30*20/(C1)*1000;
%% The main loop running the algorithm
for i = 1:T

    [d_pred, cp, sb, PV_pred] = fhandle(i, N);
    
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
            yref = [18;18;18]; %Celsius degrees -> Reference to MIN
            xref(:,i+1) = C\yref;
            ut(:,i) = Bu\(xref(:,i+1) - A*xref(:,i) - Bd*d_pred(:,1));
            %Limits on the input
            if (ut(:,i) < 0)
                ut(:,i) = 0;
            elseif(ut(:,i) > 15)
                ut(:,i) = 15;
            end
            
            %EWH temperature reference
            Tref(:,i+1) = 18; % Celsius degrees -> Reference to MIN
            Tref(:,i+1) = CO*(CO1*Tref(:,i) + CO2*d_pred(1,1) - CO3 + CO4*uet(:,i));%EWH model
            uet(:,i) = (Tref(:,i+1)/CO - CO1*Tref(:,i) - CO2*d_pred(1,1) - CO3)/CO4;
           %Limits on the input
            if (uet(:,i) < 0)
                uet(:,i) = 0;
            elseif(uet(:,i) > 4.5)
                uet(:,i) = 4.5;
            end
            
%------ %Time between 04:00 - 06:00 ---------------------------------------
        case 2
            
            %Building temperature reference
            yref = [22;22;22]; %Celsius degrees -> Reference to MAX
            xref(:,i+1) = C\yref;
            ut(:,i) = Bu\(xref(:,i+1) - A*xref(:,i) - Bd*d_pred(:,1));
            %Limits on the input
            if (ut(:,i) < 0)
                ut(:,i) = 0;
            elseif(ut(:,i) > 15)
                ut(:,i) = 15;
            end
            
            %EWH temperature reference
            Tref(:,i+1) = 22; % Celsius degrees -> Reference to MAX
            Tref(:,i+1) = CO*(CO1*Tref(:,i) + CO2*d_pred(1,1) - CO3 + CO4*uet(:,i));%EWH model
            uet(:,i) = (Tref(:,i+1)/CO - CO1*Tref(:,i) - CO2*d_pred(1,1) - CO3)/CO4;
            %Limits on the input
            if (uet(:,i) < 0)
                uet(:,i) = 0;
            elseif(uet(:,i) > 4.5)
                uet(:,i) = 4.5;
            end
            
%------ %Time between 06:00 - 18:00 ---------------------------------------
        case 3
            
            %Check the PV production
            if (PV_power > 0)
                ut(:,i) = [PV_power(i)/3;PV_power(i)/3;PV_power(i)/3];
                if ( yt(:,i) == [22;22;22])
                    
                yref = [22;22;22]; %Celsius degrees
                xref(:,i+1) = C\yref;
                ut(:,i) = Bu\(xref(:,i+1) - A*xref(:,i) - Bd*d_pred(:,1));
            end
            %Building temperature reference
            yref = [22;22;22]; %Celsius degrees
            xref(:,i+1) = C\yref;
            ut(:,i) = Bu\(xref(:,i+1) - A*xref(:,i) - Bd*d_pred(:,1));
            %Limits on the input
            if (ut(:,i) < 0)
                ut(:,i) = 0;
            elseif(ut(:,i) > 15)
                ut(:,i) = 15;
            end
            
            %EWH temperature reference
            Tref(:,i+1) = 22; % Celsius degrees 
            Tref(:,i+1) = CO*(CO1*Tref(:,i) + CO2*d_pred(1,1) - CO3 + CO4*uet(:,i));%EWH model
            uet(:,i) = (Tref(:,i+1)/CO - CO1*Tref(:,i) - CO2*d_pred(1,1) - CO3)/CO4;
            %Limits on the input
            if (uet(:,i) < 0)
                uet(:,i) = 0;
            elseif(uet(:,i) > 4.5)
                uet(:,i) = 4.5;
            end
            
%------ %Time between 18:00 - 22:00 ---------------------------------------
        case 4
            
            %Building temperature reference
            yref = [18;18;18]; %Celsius degrees
            xref(:,i+1) = C\yref;
            ut(:,i) = Bu\(xref(:,i+1) - A*xref(:,i) - Bd*d_pred(:,1));
            %Limits on the input
            if (ut(:,i) < 0)
                ut(:,i) = 0;
            elseif(ut(:,i) > 15)
                ut(:,i) = 15;
            end
            
            %EWH temperature reference
            Tref(:,i+1) = 18; % Celsius degrees 
            Tref(:,i+1) = CO*(CO1*Tref(:,i) + CO2*d_pred(1,1) - CO3 + CO4*uet(:,i));%EWH model
            uet(:,i) = (Tref(:,i+1)/CO - CO1*Tref(:,i) - CO2*d_pred(1,1) - CO3)/CO4;
            %Limits on the input
            if (uet(:,i) < 0)
                uet(:,i) = 0;
            elseif(uet(:,i) > 4.5)
                uet(:,i) = 4.5;
            end
            
    end
    
    %Logging the data
    uet(:,i) = U(end,1);
    et(:,i)  = U(end-1,1);
    vt(:,i)  = U(end-2,1);
    xbt(:,i) = xb;
    cpt(:,i) = cp(1,1);
    sbt(:,i) = sb(1,1);
    yt(:,i)  = C*xt(:,i);
    t(1,i)   = i;
    Tempt(:,i) = Temp1;
    
    disp(['Iteration ' int2str(i)]);
    yalmiperror(id);
    xt(:,i+1) = A*xt(:,i) + Bu*ut(:,i) + Bd*d_pred(:,1); %Building model equation
    if(vt(:,i)>0)
    xb = a*xb + b_ch*vt(:,i)/3; %Battery model equation
    else
    xb = a*xb + 1/b_dch*vt(:,i)/3; %Battery model equation
    end
    
    %EWH model equation
    Tempt(:,i+1) = CO*(CO1*Tempt(:,i) + CO2*d_pred(1,1) - CO3 + CO4*uet(:,i));%EWH model

end