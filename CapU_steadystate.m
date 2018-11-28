function [ys,check,xfs,data_SS,flaga]=CapU_steadystate(ys,exe)
%global alfa b betta cap chi cpp deltai epsilon eta_i g_ss gamma_p gamma_w gammaP gamma_comp h i_bar K_ss Iss kappa_pie kappa_y lambda mu_mark omega psi_ss rho rho_a rho_ksi rho_g rho_mark rho_w tau theta thetaw varphi;
check=0;
load params
options     =   optimset('TolFun',1e-28,'MaxFunEvals',2000,'Display','Off','Algorithm','trust-region-dogleg');

L_ss=0.6;
%K_ss=13;
XX=zeros(2,1);
XX(1)=K_ss;
XX(2)=L_ss;

% Calculate K_ss
[xfs, fsfv, flaga]  =   fsolve('solveSteadyStateCapU',XX,options,alfa,b,betta,cap,chi,cpp,deltai,epsilon,eta_i,g_ss,gamma_p,gamma_w,gammaP,gamma_comp,h,i_bar,K_ss,Iss,kappa_pie,kappa_y,kappa_pieW,kappa_prem,lambda,mu_mark,omega,Premium_SS,psi_ss,rho,rho_a,rho_ksi,rho_g,rho_mark,rho_w,tau,theta,thetaw,varphi);%,alfa,betta,chi,cpp,delta,epsilon,eta_i,g_ss,h,i_bar,Iss,kappa_pie,kappa_y,lambda,mu_mark,omega,psi_ss,rho,rho_a,rho_ksi,rho_g,rho_mark,rho_w,tau,theta,thetaw,varphi);

%[xfs, fsfv, flaga]  =   fsolve('solveSteadyState',XX,options,alfa,betta,chi,cpp,delta,epsilon,eta_i,g_ss,h,i_bar,Iss,kappa_pie,kappa_y,lambda,mu_mark,omega,psi_ss,rho,rho_a,rho_ksi,rho_g,rho_mark,rho_w,tau,theta,thetaw,varphi);

%% Calculate Steady State with K from fsolve
K=real(xfs(1));
L=real(xfs(2));
S=K;
Sp=S;
R=1/betta;
i=R;
Lambda=1;
A=0;
markup_w=0;
g=g_ss;
psi=0;
delta= deltai;
zeta=0;
In=0;
pieComp=1;
Q=1;
pie=1;
Pie_W=1;
pop=1;
pop_w=1;
Pm=(epsilon-1)/epsilon;
Y=K^alfa*L^(1-alfa);
G=g_ss*Y;
I=deltai*K;
C=Y-I-G;
Rk=(epsilon-1)/(epsilon)*alfa*Y/K+1-deltai;
Prem=Rk-1/betta;
varrho=(1-betta*h)/((1-h)*C);
% L=(Y/K^alfa)^(1/(1-alfa));
w=(thetaw/(thetaw-1))*chi*L^varphi/varrho;
x1_w=varrho*w*L;
x2_w=(thetaw-1)/thetaw;
x1=Y/(1-gamma_p*betta);
x2=Pm*Y/(1-gamma_p*betta);
aa     =   lambda*betta*theta*Prem;
bb     =   -(1-theta)*(lambda-betta*Prem);
cc     =   (1-theta);
phi    =   real((-bb-sqrt(bb^2-4*aa*cc))/(2*aa));
phic=phi;
psi=psi_ss;
z=(Rk-R)*phi+R;
x=z;
N=omega*K/(1-theta*(Prem*phi+1/betta));
Nn=omega*K;
Ne=N-Nn;
eta=(1-theta)/(1-betta*theta*z);
nu=((1-theta)*betta*(Rk-R))/(1-betta*theta*x);
Y_log =log(Y);
C_log =log(C);
L_log =log(L);
w_log =log(w);
Q_log =log(Q);
K_log =log(K);
I_log =log(I);
N_log =log(N);
phic_log=log(phic);
MC_log=log(Pm);
MCInv=1/Pm;
U=1;
b       =    Pm*alfa*Y/K;
delta_c=deltai-b/(1+cap);
xfs=[K,b,delta_c,Prem,L];

data_SS = [ C/Y
    I/Y
    L/Y
    w*L/Y
    %N/Y
    Prem
    Sp/Y
    z
    x];

%% Save Output
ys=[K
    S
    Sp
    R
    i
    Lambda
    A
    markup_w
    g
    psi
    zeta
    In
    pieComp
    Q
    pie
    Pie_W
    pop
    pop_w
    Pm
    Prem
    Rk
    Y
    G
    C
    varrho
    L
    I
    w
    x1_w
    x2_w
    x1
    x2
    phi
    phic
    psi
    z
    x
    N
    Nn
    Ne
    eta
    nu
    Y_log
    C_log
    L_log
    w_log
    N_log
    Q_log
    K_log
    I_log
    phic_log
    MC_log
    MCInv
    U
    b
    delta_c
    ];

%save steadyState.mat ys
end