% Thrust not taken into account, Cl of the aileron assumed from the moment
% equalibrium without thrust

Sw = 73.73;
St = 18.43;
bw = 20.95;
bt = 9.31;
cw = 3.79;
ct = 2.04;
ARw = 5.93;
ARt = 4.74;
TRw = 0.377;
TRt = 0.528;
dw = 2*pi/180;
dt = 0;
gtw = -3*pi/180;
gtt = 0;
atw = -3.5*pi/180;
att = 0;
We = 13434*0.45359237*9.81;
Wp = 1.5029e+03*0.45359237*9.81;
Wf = 5500*0.45359237*9.81;
W = (We+Wp+Wf);
alt = 3000; % meter
U = convvel(200,'kts','m/s');
[~, ~, P, rho] = atmosisa(alt);
T = 18454;
Vh = 0.6722;
lt = Vh*Sw*cw/St;
lt_dash = lt;
ht = 3.5;
% aw = 0.1154;
% at = 0.1100;
% a = aw+at*St/Sw;
a = 5.7422;
e = 1.0124;

Cl = W/(0.5*rho*Sw*U^2);
Cd = GetCd(Cl);
[eta_a_H,x_np,Cm0] = GetParYear1();
at = eta_a_H*180/pi;

x_cg_zf = 10.14;
x_cg_f = 10.06;
x_cg = (x_cg_zf*(Wp+We)+x_cg_f*Wf)/W;
Kn = (x_np-x_cg)/cw;
Clt = (Cm0-Kn*Cl)/Vh;
Cdt = 0.00540; % for Re = 10^7
d_eta_d_alfa = 1;

X_u = -rho*U*Sw*Cd;
Z_u = -rho*U*Sw*Cl;
M_u = 0;

X_w = 0.5*rho*U*Sw*(-2*a*Cl/(pi*ARw*e) + Cl);
Z_w = -0.5*rho*U*Sw*(a + Cd);
M_w = -0.5*rho*U*Sw*cw*Kn*a;

X_w_h = 0.5*rho*U*St*(-2*eta_a_H*Clt/(pi*ARt*1) + Clt);
Z_w_h = -0.5*rho*U*St*(eta_a_H + Cdt);

X_q = 0.5*rho*U*St*lt_dash*(Clt - 2*eta_a_H*Clt/(pi*ARt*1));
Z_q = -0.5*rho*U*St*lt_dash*(eta_a_H + Cdt);
M_q = lt_dash^2*Z_w_h; - ht*lt_dash*X_w_h;

X_w_dot = lt*d_eta_d_alfa*X_w_h/U;
Z_w_dot = lt*d_eta_d_alfa*Z_w_h/U;
M_w_dot = lt*d_eta_d_alfa*(-ht*X_w_h + lt_dash*Z_w_h)/U;

m0 = (We+Wp)/9.81;
mf = Wf/9.81;
m = m0+mf;
Ixx0 = m0*2.52^2;
Iyy0 = m0*4.16^2;
Izz0 = m0*4.65^2;
Ixz0 = 0.0177*m0;
Ixxf = m0*0.61^2;
Izzf = m0*0.61^2;

Ixx = Ixx0+Ixxf;
Izz = Izz0+Izzf+mf*(x_cg_f-x_cg)^2+m0*(x_cg_zf-x_cg)^2;
Iyy = Iyy0+mf*(x_cg_f-x_cg)^2+m0*(x_cg_zf-x_cg)^2;
Ixz = Ixz0;

Ms= [   m,      -X_w_dot,   0,      0;
        0,      m-Z_w_dot,  0,      0;
        0,      -M_w_dot,   Iyy,    0;
        0       0,          0,      1];

Fs = [  X_u,    X_w,    X_q,    -m*9.81;
        Z_u,    Z_w,    Z_q+m*U,0;
        M_u,    M_w,    M_q,    0;
        0,      0,      1,      0];

Y_v = 0.5*rho*U*Sw * -0.5772;
L_v = 0.5*rho*U*Sw * -0.0723;
N_v = 0.5*rho*U*Sw * 0.0773;
Y_p = 0.5*rho*U*Sw * 0.0170;
L_p = 0.5*rho*U*Sw * -0.2228;
N_p = 0.5*rho*U*Sw * -0.0121;
Y_r = 0.5*rho*U*Sw * 0.1312;
L_r = 0.5*rho*U*Sw * 0.0782;
N_r = 0.5*rho*U*Sw * -0.0645;
Y_dA = 0.5*rho*U*Sw * 0.0000;
L_dA = 0.5*rho*U*Sw * -0.0763;
N_dA = 0.5*rho*U*Sw * 0.0000;
Y_dR = 0.5*rho*U*Sw * 0.3090;
L_dR = 0.5*rho*U*Sw * 0.0250;
N_dR = 0.5*rho*U*Sw * -0.1117;

Ma = [  m,      0,      0,      0,      0;
        0,      Ixx,    -Ixz,   0,      0;
        0,      -Ixz,   Izz,    0,      0;
        0,      0,      0,      1,      0;
        0,      0,      0,      0,      1;];

Fa = [  Y_v,    Y_p,    Y_r-m*U,    m*9.81, 0;
        L_v,    L_p,    L_r,        0.      0;
        N_v,    N_p,    N_r,        0,      0;
        0,      1,      0,          0,      0;
        0,      0,      1,          0,      0];



