function [ys,check]=GK_Nom_CapU_steadystate(ys,exe)
global M_ 
%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for ii = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(ii,:));                   %    Get the name of parameter i. 
  eval([ paramname ' = M_.params(' int2str(ii) ');']);         %    Get the value of parameter i.
  eval([ 'params.' paramname '= M_.params(' int2str(ii) ');']);
end                                                           % End of the loop.  
check = 0;
%%
%% END OF THE FIRST MODEL INDEPENDENT BLOCK.
%options     =   optimset('TolFun',1e-28,'MaxFunEvals',2000,'Display','Off','Algorithm','trust-region-dogleg');

L_ss=0.5;
%K_ss=13;
XX=zeros(2,1);
XX(1)=K_ss;
XX(2)=L_ss;

% Calculate K_ss
options = optimset('TolFun',1e-28,'MaxFunEvals',2000,'Display','Off');
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
markup_p    = 0;
EWealth =0;
Eint    =0;
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
G=g*Y;
I=delta*K;
C=Y-I-G;%-G
Rk=(epsilon-1)/(epsilon)*alfa*Y/K+1-delta;
Premium=Rk-1/betta;
varrho=(1-betta*h)/((1-h)*C);
% L=(Y/K^alfa)^(1/(1-alfa));
w=Pm*(1-alfa)*(Y/L);
x1_w=varrho*w*L/(1-betta*gamma_w);
%x2_w=(thetaw-1)/thetaw*x1_w;
x2_w=chi*L^(1+varphi)/(1-betta*gamma_w);
x1=Y/(1-gamma_p*betta);
x2=Pm*Y/(1-gamma_p*betta);
aa     =   lambda*betta*theta*Premium;
bb     =   -(1-theta)*(lambda-betta*Premium);
cc     =   (1-theta);
phi    =   real((-bb-sqrt(bb^2-4*aa*cc))/(2*aa));
phic=phi;
psi=psi_ss;
z=(Rk-R)*phi+R;
x=z;
N=omega*K/(1-theta*(Premium*phi+1/betta));
Nn=omega*K;
Ne=N-Nn;
eta=(1-theta)/(1-betta*theta*z);
nu=((1-theta)*betta*(Rk-R))/(1-betta*theta*x);
MPK = alfa*Y/K;
MPL = (1-alfa)*Y/L;
MPK_log = 0;
MPL_log = 0;
Y_log =0; 
C_log =0;	
L_log =0;
w_log =0; 
Q_log =0; 
K_log =0; 
I_log =0; 
N_log =0;
phic_log=0;
MC_log=0; 
MCInv=1/Pm;
U=1;
b       =    Pm*alfa*Y/K;  
delta_c=delta-b/(1+cap);
Incentive = (nu*(Q*S)+eta*N) - lambda*(Q*S);

eps_a=0;		
ksi=0;
eps_m =0;
eps_w=0;
eps_g=0;
eps_mark =0;

varrho_Flex = varrho;
C_Flex 		= C;
Lambda_Flex = Lambda;
R_Flex 		= R;
w_Flex 		= w;
L_Flex 		= L;
Y_Flex 		= Y;
I_Flex 		= I;
delta_Flex  = delta;
K_Flex 		= K;
In_Flex 	= In;
psi_Flex 	= psi;
Q_Flex 		= Q;
S_Flex 		= S;
Premium_Flex= Premium; 
phi_Flex  	= phi;
phic_Flex 	= phic;
N_Flex 		= N;
Sp_Flex 	= Sp;
U_Flex 		= U;
z_Flex 		= z;
x_Flex 		= x;
Rk_Flex 	= Rk;
eta_Flex 	= eta;
nu_Flex     = nu;
G_Flex      = G;
y_gap 		= Y/Y_Flex;
%y_gap2 		= (epsilon/(epsilon-1))/MCInv;
y_gap2      = Y-Y_Flex;
Wf 			= (log(C*(1-h))-chi/(1+varphi)*L^(1+varphi))/(1-betta); 
EUL 		= chi/(1+varphi)*L^(1+varphi)/(1-betta); 
C_lag		= C;
Rk_nom      = Rk;
x_obs       = 0;
Wf_obs      = 0;
EPi         = 0;
%Wf_obs      = log(Wf);
%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state values of the endogenous variables of
%% the model.
%%
NumberOfEndogenousVariables = M_.orig_endo_nbr;                % Number of endogenous variables.
ys = zeros(NumberOfEndogenousVariables,1);                     % Initialization of ys (steady state).
for ii = 1:NumberOfEndogenousVariables                         % Loop...
  varname = deblank(M_.endo_names(ii,:));                      %    Get the name of endogenous variable i.                     
  eval(['ys(' int2str(ii) ') = ' varname ';']);                %    Get the steady state value of this variable.
end
end
