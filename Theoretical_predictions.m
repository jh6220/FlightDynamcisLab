% should IAS or TAS be used?
g = 9.81;
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
U = convvel(223,'kts','m/s');
[~, ~, P, rho] = atmosisa(alt);
T = 18454;
Vh = 0.6722;
lt = Vh*Sw*cw/St;
lt_dash = lt;
ht = 3.5;
lT = 0.28; % Thurst line offset
% aw = 0.1154;
% at = 0.1100;
% a = aw+at*St/Sw;
a = 5.7422;
e = 1.0124;
alpha0 = 0;

Cl = W/(0.5*rho*Sw*U^2);
Cd = GetCd(Cl);
[eta_a_H,x_np,Cm0] = GetParYear1();
at = eta_a_H*180/pi;

x_cg_zf = 10.14;
x_cg_f = 10.06;
x_cg = (x_cg_zf*(Wp+We)+x_cg_f*Wf)/W;
Kn = (x_np-x_cg)/cw;
Clt = (Cm0-Kn*Cl)/Vh;
Cdt0 = 0.00540; % for Re = 10^7
Cdt = Cdt0 + Clt^2/(pi*ARt*1);
d_eta_d_alfa = 1;

X_u = -rho*U*Sw*Cd;
Z_u = -rho*U*Sw*Cl;
M_u = 0;

X_w = 0.5*rho*U*Sw*(-2*a*Cl/(pi*ARw*e) + Cl);
Z_w = -0.5*rho*U*Sw*(a + Cd);
M_w = -0.5*rho*U*Sw*cw*Kn*a;

X_w_h = 0.5*rho*U*St*(-2*at*Clt/(pi*ARt*1) + Clt);
Z_w_h = -0.5*rho*U*St*(at + Cdt);

X_q = 0.5*rho*U*St*lt_dash*(Clt - 2*at*Clt/(pi*ARt*1));
Z_q = -0.5*rho*U*St*lt_dash*(at + Cdt);
M_q = lt_dash^2*Z_w_h - ht*lt_dash*X_w_h;

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
L_v = 0.5*rho*U*Sw*bw * -0.0723;
N_v = 0.5*rho*U*Sw*bw * 0.0773;
Y_p = 0.5*rho*U*Sw*bw * 0.0170;
L_p = 0.5*rho*U*Sw*bw^2 * -0.2228;
N_p = 0.5*rho*U*Sw*bw^2 * -0.0121;
Y_r = 0.5*rho*U*Sw*bw * 0.1312;
L_r = 0.5*rho*U*Sw*bw^2 * 0.0782;
N_r = 0.5*rho*U*Sw*bw^2 * -0.0645;
Y_dA = 0.5*rho*U*Sw * 0.0000;
L_dA = 0.5*rho*U*Sw*bw * -0.0763;
N_dA = 0.5*rho*U*Sw*bw * 0.0000;
Y_dR = 0.5*rho*U*Sw * 0.3090;
L_dR = 0.5*rho*U*Sw*bw * 0.0250;
N_dR = 0.5*rho*U*Sw*bw * -0.1117;

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

[vecs,vals] = eig(Ms\Fs);
% for phugoid predicts non-oscilatory behaviour

% disp(vals(1,1))

% SPPO
Msppo = [   m-Z_w_dot,     0;
            -M_w_dot,       Iyy];
Fsppo = [   Z_w,        Z_q+m*U;
            M_w,        M_q];

[vecsppo,valsppo] = eig(Msppo\Fsppo);

% Phugoid
Mph = [m,      -X_w_dot;
       0,      m-Z_w_dot];
Fph = [X_u,    X_w;
       Z_u,    Z_w];

[vecph,valph] = eig(Mph\Fph); % for phugoid predicts non-oscilatory behaviour

gamma = -Z_w/m;
delta = -X_w*Z_u/m^2;
lambdaph = -gamma/2 + ([1,-1]*sqrt((gamma/2)^2 - 2*g^2/U^2)); % for phugoid predicts non-oscilatory behaviour

% Assymetric
[veca,vala] = eig(Ma\Fa);

% Roll substidence
lambda_roll = L_p/Ixx;

% Spiral mode
lambda_spiral = -g*(L_v*N_r-L_r*N_v)/(U*(L_v*N_p-L_p*N_v));

% Dutch roll
lambda_dr = 0.5*(N_r/Izz+Y_v/m) + [1,-1]*sqrt(0.25*(N_r/Izz+Y_v/m)^2 - (Y_v*N_r-Y_r*N_v)/(m*Izz) - (U*N_v)/Izz);

lambda_dr2 = 0.5*(N_r/Izz+Y_v/m) + [1,-1]*sqrt(-U*N_v/Izz);

X_dE = -0.5*rho*U^2*St*(2*Clt*at/(pi*ARt*1));
Z_dE = -0.5*rho*U^2*St*at;
M_dE = lt_dash*Z_dE - ht*X_dE;

X_T = cos(alpha0);
Z_T = -sin(alpha0);
M_T = lT;

% Longitudinal data for export
long.m = m;
long.Iyy = Iyy;
long.Ue = U;

long.Xu = X_u;
long.Zu = Z_u;
long.Mu = M_u;

long.Xw = X_w;
long.Zw = Z_w;
long.Mw = M_w;

long.Xq = X_q;
long.Zq = Z_q;
long.Mq = M_q;

long.Xwdot = X_w_dot;
long.Zwdot = Z_w_dot;
long.Mwdot = M_w_dot;

long.XiH = X_dE;
long.ZiH = Z_dE;
long.MiH = M_dE;

long.Xt = X_T;
long.Zt = Z_T;
long.Mt = M_T;