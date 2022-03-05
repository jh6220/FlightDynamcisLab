function [eta_a_H,x_np,Cm0] = GetParYear1()
alt = distdim(5000,'ft','m');
[~,~,~,rho] = atmosisa(alt);
[~,~,~,rho0] = atmosisa(0);
S_ref = 73.73;
V_H = 0.6305;
cord = 3.79;
S_H = 18.43;
cg_ZF = 10.14;
cg_F = 10.06;
EW = convforce(13434.8,'lbf','N');
PW = convforce(2003.92,'lbf','N');

U_i_H = 0.05;
[~,~,~,rho1] = atmosisa(distdim(4500,'ft','m'));
u_rho = abs(rho-rho1)/rho;
U_IAS = 1;
U_W = 50;

IAS0 = convvel([250,210,170,130,90,80],'kts','m/s');
FW0 = convforce([6430,6408,6395,6379,6363,6350],'lbf','N');
AoA0 = [-1	-0.4	0.65	2.75	8	10.2];
i_H0 = [0.74	0.355	0.035	-0.915	-2.7	-3.3];
Cl0 = (EW+PW+FW0)./(0.5.*rho0.*IAS0.^2.*S_ref);
off0 = zeros(1,length(IAS0));
lin0 = polyfit(Cl0,i_H0,1);
m0 = fitlm(Cl0,i_H0);
r0 = m0.Rsquared.Ordinary.*100;
U_Cl0 = Cl0.*sqrt(u_rho.^2+2.*(U_IAS./IAS0).^2+(U_W./(EW+PW+FW0)).^2);

IASp = convvel([249,210,170,130,90,80],'kts','m/s');
FWp = convforce([6253,6236,6223,6210,6194,6183],'lbf','N');
AoAp = [-1	-0.4	0.6	2.7	7.7	10.1];
i_Hp = [0.83	0.53	0.1	-0.575	-1.95	-2.62];
Clp = (EW+PW+FWp)./(0.5.*rho0.*IASp.^2.*S_ref);
offp = convlength(6.4,'in','m').*ones(1,length(IASp));
linp = polyfit(Clp,i_Hp,1);
mp = fitlm(Clp,i_Hp);
rp = mp.Rsquared.Ordinary.*100;
U_Clp = Clp.*sqrt(u_rho.^2+2.*(U_IAS./IASp).^2+(U_W./(EW+PW+FWp)).^2);

IASn = convvel([250,211,170,130,92,81],'kts','m/s');
FWn = convforce([6015,6000,5980,5970,5944,5928],'lbf','N');
AoAn = [-1	-0.4	0.65	2.8	8	10.1];
i_Hn = [0.62	0.25	-0.275	-1.165	-3.12	-3.9];
Cln = (EW+PW+FWn)./(0.5.*rho0.*IASn.^2.*S_ref);
offn = -convlength(6,'in','m').*ones(1,length(IASn));
linn = polyfit(Cln,i_Hn,1);
mn = fitlm(Cln,i_Hn);
rn = mn.Rsquared.Ordinary.*100;
U_Cln = Cln.*sqrt(u_rho.^2+2.*(U_IAS./IASn).^2+(U_W./(EW+PW+FWn)).^2);

% errorbar(Cl0,i_H0,U_i_H*ones(1,length(i_H0)),U_i_H*ones(1,length(i_H0)),U_Cl0,U_Cl0,'xr')
% hold on
% plot(Cl0,lin0(2)+lin0(1).*Cl0,'--r')
% hold on
% errorbar(Clp,i_Hp,U_i_H*ones(1,length(i_Hp)),U_i_H*ones(1,length(i_Hp)),U_Clp,U_Clp,'xb')
% hold on
% plot(Clp,linp(2)+linp(1).*Clp,'--b')
% hold on
% errorbar(Cln,i_Hn,U_i_H*ones(1,length(i_Hn)),U_i_H*ones(1,length(i_Hn)),U_Cln,U_Cln,'xg')
% hold on
% plot(Cln,linn(2)+linn(1).*Cln,'--g')
% legend("0 in offset","0 in offset linear fit (R^2 = "+num2str(r0)+"%)","+6.4 in offset","+6.4 in offset linear fit (R^2 = "+num2str(rp)+"%)", "-6 in offset", "-6 in offset linear fit (R^2 = "+num2str(rn)+"%)")
% xlabel("Cl")
% ylabel("i_H (deg)")

Cl = [Cl0,Clp,Cln];
i_H = [i_H0,i_Hp,i_Hn];
FW = [FW0,FWp,FWn];
offset = [off0,offp,offn];
x_cg = ((EW+PW).*(cg_ZF + offset) + FW.*cg_F)./(EW+PW+FW);

A = [1./V_H.*ones(1,length(Cl));Cl./(V_H.*cord);(x_cg.*Cl)./(V_H.*cord)];
A = transpose(A);
i_H = transpose(i_H);

b = regress(i_H,A);

x1 = linspace(min(A(:,2)),max(A(:,2)),100);
x2 = linspace(min(A(:,3)),max(A(:,3)),100);
[X1,X2] = meshgrid(x1,x2);

Y = b(1) + b(2).*X1 + b(3).*X2;

% scatter3(A(:,2),A(:,3),i_H,'filled')
% hold on
% mesh(X1,X2,Y)
% zlim([-4,2])
% xlabel('X1')
% ylabel('X2')
% zlabel('Y')

eta_a_H = 1/b(3);

x_np = -b(2)*eta_a_H;

Cm0 = b(1).*eta_a_H;
end