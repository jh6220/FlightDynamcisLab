%Code to read the output data file from the Simulink recorder
%Imperial College London
%Flight Sim Lab 20/21

clear
clc
close all

%%%%%INSERT THE FILENAME HERE%%%%%
fname = 'roll_test_velocitytest.mat';
load(fname) %Load the file name
% data = sppo3;

% data.Data = Datas;
%Extract all the data
vInd_kias   = data.Data(1:100:end,1); %Indicated Airspeed - Knots
vTrue_ktas  = data.Data(1:100:end,2); %true Airspeed      - Knots
climb_rate  = data.Data(1:100:end,3); %climb rate         - fpm
q           = data.Data(1:100:end,7); %Pitch rate         - rad/s  
p           = data.Data(1:100:end,8); %Roll rate          - rad/s
r           = data.Data(1:100:end,9); %Yaw rate           - rad/s
pitch       = data.Data(1:100:end,10);%Pitch angle        - deg
roll        = data.Data(1:100:end,11);%Roll angle         - deg
heading_true= data.Data(1:100:end,12);%Heading            - deg
alpha       = data.Data(1:100:end,13);%Angle of Attack    - deg
beta        =data.Data(1:100:end,14); %Side slip angle    - deg
latitude    = data.Data(1:100:end,15);%Latitude           - deg
longitude   = data.Data(1:100:end,16);%Longitude          - deg  
altitude    = data.Data(1:100:end,17);%Altitude           - ft 
x           = data.Data(1:100:end,18);
y           = data.Data(1:100:end,19);
z           = data.Data(1:100:end,20);
throttle_cmd = data.Data(1:100:end,21);%Throttle command  - %
throttle_actual = data.Data(1:100:end,22);
eng_power   = data.Data(1:100:end,23); %Engine power      - hp
w_empty     = data.Data(1:100:end,24); 
w_payld     = data.Data(1:100:end,25); %Payload weight    - lb
w_fuel      = data.Data(1:100:end,26); %Fuel weight       - lb 
time        = data.Data(1:100:end,27); %Time from start   - seconds
elevator    = data.Data(1:100:end,28); %Tail incidence    - % 
aileron     = data.Data(1:100:end,29); %                  - %
rudder      = data.Data(1:100:end,30); %                  - %
t = data.Time(1:100:end,1);