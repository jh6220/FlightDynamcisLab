alt = distdim(5000,'ft','m');
[T,a,P,rho] = atmosisa(alt);
[T0,a0,P0,rho0] = atmosisa(0);
S_ref = 73.73;
EW = convforce(13434.8,'lbf','N');
PW = convforce(2003.92,'lbf','N');

FWs_clean = convforce(3900,'lbf','N');
vs_clean = convvel(57,'kts','m/s');
vsH_clean = convvel(mean([75,78]),'kts','m/s');

FWs_landing = convforce(3755,'lbf','N');
vs_landing = convvel(51,'kts','m/s');
vsH_landing = convvel(61,'kts','m/s');

Cl_max_clean = 2*(EW+PW+FWs_clean)/(rho0*vs_clean^2*S_ref)
Cl_max_landing = 2*(EW+PW+FWs_landing)/(rho0*vs_landing^2*S_ref)

glide_T = [12.59,15.8,20.3,26.3,38.6,48.7,62.7,75.7,82.6,94.5,92.8];
glide_IAS = convvel([287,266,245,225,200,180,160,145,130,115,95],'kts','m/s');
glide_FW = convforce([5910,5860,5810,5762,5928,5888,5845,5790,5713,5661,5911],'lbf','N');

Cl = 2.*(EW+PW+glide_FW)./(rho0.*glide_IAS.^2.*S_ref);
Cd = Cl.*sqrt(rho/rho0).*distdim(1000,'ft','m')./(glide_IAS.*glide_T);

scatter(Cd,Cl)
xlabel("Cd")
ylabel("Cl")
p = polyfit(Cl,Cd,2)
P = polyfit(Cl.^2,Cd,1)
Y = linspace(0,max(Cl),100);
X = p(1).*Y.^2+p(2).*Y+p(3);
hold on
plot(X,Y,'--')
X = P(1).*Y.^2+P(2);
% P(1) = 1/(AR*pi*e)
AR = 5.93;
e = 1/(pi*AR*P(1))
hold on
plot(X,Y)
legend("experimental data", "quadratic best fit")



