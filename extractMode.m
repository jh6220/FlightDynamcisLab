function [mode] = extractMode(fname)

data = load(fname);

mode.vInd_kias   = data.Datas(:,1); %Indicated Airspeed - Knots
mode.vTrue_ktas  = data.Datas(:,2); %true Airspeed      - Knots
mode.climb_rate  = data.Datas(:,3); %climb rate         - fpm
mode.q           = data.Datas(:,4); %Pitch rate         - rad/s  
mode.p           = data.Datas(:,5); %Roll rate          - rad/s
mode.r           = data.Datas(:,6); %Yaw rate           - rad/s
mode.pitch       = data.Datas(:,7);%Pitch angle        - deg
mode.roll        = data.Datas(:,8);%Roll angle         - deg
mode.heading_true= data.Datas(:,9);%Heading            - deg
mode.alpha       = data.Datas(:,10);%Angle of Attack    - deg
mode.beta        =data.Datas(:,11); %Side slip angle    - deg
mode.latitude    = data.Datas(:,12);%Latitude           - deg
mode.longitude   = data.Datas(:,13);%Longitude          - deg  
mode.altitude    = data.Datas(:,14);%Altitude           - ft 
mode.x           = data.Datas(:,15);
mode.y           = data.Datas(:,16);
mode.z           = data.Datas(:,17);
mode.throttle_cmd = data.Datas(:,18);%Throttle command  - %
mode.throttle_actual = data.Datas(:,19);
mode.eng_power   = data.Datas(:,20); %Engine power      - hp
mode.w_empty     = data.Datas(:,21); 
mode.w_payld     = data.Datas(:,22); %Payload weight    - lb
mode.w_fuel      = data.Datas(:,23); %Fuel weight       - lb 
mode.time        = data.Datas(:,24); %Time from start   - seconds
mode.elevator    = data.Datas(:,25); %Tail incidence    - % 
mode.aileron     = data.Datas(:,26); %                  - %
mode.rudder      = data.Datas(:,27); %                  - %

end