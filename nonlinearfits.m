%nlinfits
%%

clear
clc
close all

fname = ['LeftRoll1.mat'];
data = load(fname);

vInd_kias   = data.Datas(:,1); %Indicated Airspeed - Knots
vTrue_ktas  = data.Datas(:,2); %true Airspeed      - Knots
climb_rate  = data.Datas(:,3); %climb rate         - fpm
q           = data.Datas(:,4); %Pitch rate         - rad/s  
p           = data.Datas(:,5); %Roll rate          - rad/s
r           = data.Datas(:,6); %Yaw rate           - rad/s
pitch       = data.Datas(:,7);%Pitch angle        - deg
roll        = data.Datas(:,8);%Roll angle         - deg
heading_true= data.Datas(:,9);%Heading            - deg
alpha       = data.Datas(:,10);%Angle of Attack    - deg
beta        =data.Datas(:,11); %Side slip angle    - deg
latitude    = data.Datas(:,12);%Latitude           - deg
longitude   = data.Datas(:,13);%Longitude          - deg  
altitude    = data.Datas(:,14);%Altitude           - ft 
x           = data.Datas(:,15);
y           = data.Datas(:,16);
z           = data.Datas(:,17);
throttle_cmd = data.Datas(:,18);%Throttle command  - %
throttle_actual = data.Datas(:,19);
eng_power   = data.Datas(:,20); %Engine power      - hp
w_empty     = data.Datas(:,21); 
w_payld     = data.Datas(:,22); %Payload weight    - lb
w_fuel      = data.Datas(:,23); %Fuel weight       - lb 
time        = data.Datas(:,24); %Time from start   - seconds
elevator    = data.Datas(:,25); %Tail incidence    - % 
aileron     = data.Datas(:,26); %                  - %
rudder      = data.Datas(:,27); %                  - %


%%
%-------------

%%%%%%%%%%%%%%%%%%%%%%%% OSCILLATORY
%modelfun = @(k,t) k(1) + k(2).*exp(k(3).*t).*sin((k(4).*t + k(5))); %Oscillatory

%%%%%%%%%%%%%%%%%%%%%%%% NON OSCILLATORY
modelfun = @(k,t) k(1) + k(2).*exp(k(3).*t); %Non Oscillatory
figure; hold on;


Y = roll(42458:84628);
timeToUSE = time(42458:84628);
Y = Y(1:32194);
timeToUSE = timeToUSE(1:32194);
timeToUSE = timeToUSE-min(timeToUSE);
plot(timeToUSE,Y,'r');

%%

figure; hold on
% bC1 = 0.5
% bC2 = 0.5
% bC3 = -0.5
% bC4 = 0.5
% bC5 = 0.5
% betaCoefficient0 = [bC1 bC2 bC3 bC4 bC5]'; %Oscillatory

bC1 = 0.5;
bC2 = 0.5;
bC3 = 0.5;
betaCoefficient0 = [bC1 bC2 bC3]'; %NonOscillatory


betaCoefficient = nlinfit(timeToUSE,Y,modelfun,betaCoefficient0);


Ycoords = modelfun(betaCoefficient,timeToUSE);
plot(timeToUSE,Ycoords,'b');


figure; hold on
plot(timeToUSE,Y,'r');
plot(timeToUSE,Ycoords,'b');

%
%%
for i = 1:27
    
    if i ~= 24
        figure
        plot(time,data.Datas(:,i),'Color',[rand rand rand])
        switch i
            case 1
                word = 'vInd_kias';
            case 2
                word = 'vTrue_ktas';
            case 3
                word = 'climb_rate';
            case 4
                word = 'q';
            case 5
                word = 'p';
            case 6
                word = 'r';
            case 7
                word = 'pitch';
            case 8 
                word = 'roll';
            case 9
                word = 'heading_true';
            case 10
                word = 'alpha';
            case 11
                word = 'beta';
            case 12
                word = 'latitude';
            case 13
                word = 'longitude';
            case 14
                word = 'altitude';
            case 15
                word = 'x';
            case 16
                word = 'y';
            case 17
                word = 'z';
            case 18
                word = 'throttle_cmd';
            case 19
                word = 'throttle_actual';
            case 20
                word = 'eng_power';
            case 21
                word = 'w_empty';
            case 22
                word = 'w_payld';
            case 23
                word = 'w_fuel';
            case 25
                word = 'elevator';
            case 26
                word = 'aileron';
            case 27
                word = 'rudder';
        end
        title(word)
    end
end

function RandomNumber = randomInt(min, max)

    RandomNumber = min + (max-min)*rand;

end
