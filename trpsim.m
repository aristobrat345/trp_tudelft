
%Constants
g=9.81; %m/s^2
R_a=8.314; %J/mol/K 
alpha=0.861;
%Rocket parameters 
d_t=0.0112; %m +-0.00005 change?
d_e=0.022; %m +-0.00005 change?
d_out=0.036;
l_initial=0.060; %m
d_port_initial=0.020; %m
a=5.5;
n=-0.013;
T_c=1575; %K
M_w=39.5324/1000; %kg/mol
R=R_a/M_w; 
rho_s=1841;%Change this later
%S_i=%%%

%Density of the gas

%rho_g=

%Calculating exit pressure
A_e=pi*d_e^2/4;
A_t=pi*d_t^2/4;
gamma=1.1375; %This is a random value. Change this later
GAMMA=sqrt(gamma)*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)));
A_ratio=A_e/A_t; %m^2
p_ratio=0.05162; %Pa
% for p_ratio=1:10000
%     A_ratio(1,p_ratio)=GAMMA/sqrt(2*gamma/(gamma-1)*(p_ratio/10000)^(2/gamma)*(1-(p_ratio/10000)^((gamma-1)/gamma)))
%     %if abs(A_ratio-A_e/A_t)<0.01
%      %   break
%     %end
% end
% A_ratio_true=A_e/A_t*ones(1,10000);
% x=abs(A_ratio_true-A_ratio);
% y=find(min(x));

%% Initial Conditions 
c_initial=1/GAMMA*sqrt(R*T_c);
S_initial=l_initial*4*d_port_initial*pi+(d_out^2-d_port_initial^2)*0.25*pi*8;
p_c_initial=(c_initial*a/1000*rho_s*S_initial/A_t)^(1/(1-n));
rho_c=p_c_initial/R/T_c;
p_c_initial_a=(c_initial*a/1000*(rho_s-rho_c)*S_initial/A_t/alpha)^(1/(1-n));
r_initial=a/1000*p_c_initial_a^n;

%% Iteration over time
t=


%r=a*p_c^n
%d_port=d_port_initial-2*r*t;
%l_bates=l_initial-2*r*t;
%S=l_bates*4*d_port*pi+(d_out^2-d_port^2)*0.25*pi*8;
%m=rho_s*r*S
%v_e=sqrt(2*gamma/(gamma-1)*R*T_c*(1-(p_e/p_c)^((gamma-1)/gamma)))
%F_T=m*v_e+(p_e-p_a)*A_e
%c=p_c*A_t/m
%C_f=F_T/p_c/A_t
%Isp=F_T/g/m+(p_e-p_a)*A_e
