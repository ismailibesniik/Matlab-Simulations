%% R-C modeling of the Thermal Capacity of a Building
%Puprose: Master Thesis Project
%Author: Besnik Ismaili
%Script that fits the given state space model of a building to an first
%order R-C model
clc;
close all;
yalmip('clear')
clear all
load hist_weather.mat
t_20 = meteo.datetime(1):minutes(20):meteo.datetime(end);
% temp=double(str2num(char(meteo.temp)));
% 
% for i=1:length(meteo.temp)
%     temp(i) = temp(i) + rand/100;
% end
%  temp_20 = interp1(meteo.datetime,unique(temp), t_20);
%  plot(meteo.datetime(1:100),temp_20, t_20(1:50), temp_20(1:50));

d_pred = [double(str2num(char(meteo.temp)))';double(str2num(char(meteo.glob_irr)))';double(zeros(size(meteo.temp')))];