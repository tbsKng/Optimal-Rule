var
Y 			% Output
C 			% Consumption
L  			% Labour Hours
w  			% Real Wage
Pm 			% Price of intermediate Goods = MC
MCInv  		% Price Markup
R    		% Risk-free rate
i 			% Mominal Interest Rate
pie    		% CPI Inflation
pop 		% Relative Price
x1   		% Recursive Retailer Variable
x2 			% Recursive Retailer Variable 2
Pie_W  		% Wage Inflation
pop_w 		% Relative Wage
x1_w 		% Recursive Wage Variable
x2_w 		% Recursive Wage Variable 2
varrho 		% Marginal Utility
Lambda 		% Stochastic Discount Factor
A  			% Technology
I 			% Investment
In 			% Net Investment
K 			% Capital
Q       	% Price of Capital
zeta 		% Capital Quality
Rk 			% Return on Capital
nu 			% Marginal Gain by Issuing one Further Unit of Assets
eta 		% Marginal Gain by Increasing Net Wealth by one Unit
N  			% Net Wealth
Nn 			% Net Wealth of newly entering bankers
Ne 			% Net Wealth of existing banks
S 			% Assets=loans
phi 		% Leverage Ratio of private banks
z 			% Growth Rate of Net Wealth
x  			% Growth Rate of Assets
Premium 	% Premium
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
MPL_log
MPK_log
G 			% Government Expenditures
g 			% Government Expenditures relative to output
phic 		% Private+Public Leverage Ratio
psi 		% Fraction of Government Loans
Sp 			% Government Loans
markup_w 	% wage mark-up
markup_p    % price mark-up
EWealth     % Wealth shock
Eint        % Interest Rate Shock
pieComp
MPL
MPK
%MRS
delta 		% depreciation rate
U
y_gap
varrho_Flex
C_Flex
Lambda_Flex
R_Flex
w_Flex
L_Flex
Y_Flex
I_Flex
delta_Flex
K_Flex
In_Flex
psi_Flex
Q_Flex
S_Flex
Premium_Flex
phi_Flex
phic_Flex
N_Flex
Sp_Flex
U_Flex
z_Flex
x_Flex
Rk_Flex
eta_Flex
nu_Flex
G_Flex
y_gap2
Wf
C_lag
EUL
Rk_nom
x_obs
Wf_obs
Incentive
EPi
;

varexo
eps_a  		% Technology Shock
ksi 		% Capital Quality Shock
eps_m 		% Monetary Policy Shock
%eps_w 		% Wealth Shock
eps_g 		% Government expenditure shock
eps_mark_w 	% Shock to Wage Mark-up
eps_mark_p  % Shock to price mark-up
%eps_pie     % Shock to inflation rate
;


parameters
alfa 		% capital share
b
betta 		% Household's discount rate
chi 		% Disutility from Labor
cap         % Utilization Parameter
cpp 		% Credit policy parameter measuring degree of intervention
%delta
deltai
delta_c
epsilon 	% elasticity of substitution between two goods
eta_i 		% parameter for adjustment cost function
g_ss 		% Share of governemnt expenditures relative to GDP
gamma_p 		% Calvo parameter for price setting
gamma_comp
gamma_w 	% Calvo paramter for wage setting
gammaP 		% Inflation Indexation parameter
gammaW      % Wage inflation indexation
h 			% habit formation
i_bar		% St. St. nominal interest rate
Iss 		% Steady State Investment
K_ss		% St.St. capital
%L_ss
kappa_pie 	% Taylor parameter for inflation
kappa_y 	% Taylor parameter for output gap
kappa_x
kappa_prem
kappa_pieW
lambda 		% probility bankers divert assets
mu_mark 	% MA(1) coefficient for wage mark up shock
Premium_SS
%pi_ss
%piW_ss
psi_ss 		% Steady state share of government credit
omega 		% fraction of total assets transfered to newly entering bankers
rho 		% Interest Rate smooting parameter
rho_a  		% Technology autoregressive parameter
rho_Eint    % Interest rate autoregressive parameter
rho_g 		% Government expenditures autoregressive parameter
rho_ksi  	% Capital Quality autoregressive parameter
rho_mark	% AR(1) coefficient for wage mark-up shock
rho_mark_p	% AR(1) coefficient for price mark-up shock
rho_w 		% Wealth autoregressive parameter
rho_prem
sigma_epsiA
sigma_epsiK
sigma_epsiMW
sigma_epsiMP
sigma_epsiW
sigma_epsiInt
sigma_epsiPrem
%sigma_epsiPrem_obs
sigma_epsiG
tau 		% Efficiency costs of government credit policy
thetaw      % elasticity of substitution between different kind of labor
theta       % survival rate of bankers
varphi		% Inverse Frisch Labor elasticity
rho_pie
;

load params.mat
%load paramsEst.mat
set_param_value('alfa',alfa);
set_param_value('b',b);
set_param_value('betta',betta);
set_param_value('cap',cap);
set_param_value('chi',chi);
set_param_value('cpp',cpp);
set_param_value('deltai',deltai);
set_param_value('delta_c',delta_c);
set_param_value('epsilon',epsilon);
set_param_value('eta_i',eta_i);
set_param_value('g_ss',g_ss);
set_param_value('gamma_p',gamma_p);
set_param_value('gamma_w',gamma_w);
set_param_value('gammaP',gammaP);
set_param_value('gammaW',gammaW);
set_param_value('gamma_comp',gamma_comp);
set_param_value('h',h);
set_param_value('i_bar',i_bar);
set_param_value('K_ss',K_ss);
%set_param_value('L_ss',L_ss);
Iss=deltai*K_ss;
set_param_value('kappa_pie',kappa_pie);
set_param_value('kappa_y',kappa_y);
set_param_value('kappa_x',kappa_x);
set_param_value('kappa_prem',kappa_prem);
set_param_value('kappa_pieW',kappa_pieW);
set_param_value('lambda',lambda);
%set_param_value('pi_ss',pi_ss);
%set_param_value('piW_ss',piW_ss);
set_param_value('psi_ss',psi_ss);
set_param_value('Premium_SS',Premium_SS);
%set_param_value('sigma_epsiPrem_obs',sigma_epsiPrem_obs);
%set_param_value('mu_mark',mu_mark);
set_param_value('omega',omega);
psi_ss=0;
set_param_value('rho',rho);
set_param_value('rho_a',rho_a);
set_param_value('rho_Eint',rho_Eint);
set_param_value('rho_g',rho_g);
set_param_value('rho_ksi',rho_ksi);
set_param_value('rho_mark',rho_mark);
set_param_value('rho_mark_p',rho_mark_p);
set_param_value('rho_prem',rho_prem);
set_param_value('rho_w',rho_w);
set_param_value('sigma_epsiA',sigma_epsiA);
set_param_value('sigma_epsiK',sigma_epsiK);
set_param_value('sigma_epsiMW',sigma_epsiMW);
set_param_value('sigma_epsiMP',sigma_epsiMP);
set_param_value('sigma_epsiW',sigma_epsiW);
set_param_value('sigma_epsiInt',sigma_epsiInt);
set_param_value('sigma_epsiPrem',sigma_epsiPrem);
set_param_value('sigma_epsiG',sigma_epsiG);
tau=0.001;
set_param_value('theta',theta);
set_param_value('thetaw',thetaw);
set_param_value('varphi',varphi);
rho_pie = 0;
model;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Flex Price Economy %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Households
// 1
varrho_Flex=(C_Flex-h*C_Flex(-1))^(-1)-betta*h*(C_Flex(+1)-h*C_Flex)^(-1);
// 2
Lambda_Flex=varrho_Flex(+1)/varrho_Flex;
// 3
1=betta*(varrho_Flex(+1)/varrho_Flex)*R_Flex;

% Wage Rigidities
// 4
varrho_Flex*w_Flex*L_Flex = (thetaw/(thetaw-1))*exp(markup_w)*chi*L_Flex^(varphi+1);
%w_Flex= (Pie_W_Flex/pie_Flex)*w_Flex(-1);


% Economy wide resource constraint
// 5
Y_Flex=C_Flex+I_Flex+(eta_i/2)*((I_Flex-delta_Flex*exp(zeta)*K_Flex(-1)+Iss)/(In_Flex(-1)+Iss)-1)^2*(I_Flex-delta_Flex*exp(zeta)*K_Flex+Iss)+tau*psi_Flex*Q_Flex*S_Flex+G_Flex;%; %
// 6
G_Flex=g*Y_Flex;

% Credit Policy of Central Bank
// 7
psi_Flex=psi_ss+cpp*(Premium_Flex-Premium_SS);
// 8
Q_Flex*S_Flex=phic_Flex*N_Flex;
// 9
S_Flex=Sp_Flex+psi_Flex*S_Flex;
// 10
phic_Flex=phi_Flex/(1-psi_Flex);

% Intermediate Goods-producing Firms
// 11
K_Flex=S_Flex;
// 12
Y_Flex=exp(A)*(U_Flex*exp(zeta)*K_Flex(-1))^(alfa)*L_Flex^(1-alfa);   % with U
// 13
w_Flex=((epsilon-1)/epsilon)*(1-alfa)*(Y_Flex/L_Flex);							% comp labor market equilibrium
// 14
((epsilon-1)/epsilon)*exp(markup_p)=exp(zeta)*K_Flex(-1)*b*U_Flex^cap*(U_Flex/(alfa*Y_Flex));
// 15
delta_Flex = delta_c+b/(1+cap)*(U_Flex^(1+cap)); % depreciation rate

% Financial Intermediaries
// 16
N_Flex=theta*z_Flex*N_Flex(-1)*exp(EWealth)+omega*Q_Flex*S_Flex(-1);%-eps_w;
%Q*S=phi*N;
// 17
z_Flex=(Rk_Flex-R_Flex(-1))*phi_Flex(-1)+R_Flex(-1);
// 18
x_Flex=(phi_Flex/phi_Flex(-1))*z_Flex;
// 19
eta_Flex=(1-theta)+betta*Lambda_Flex*theta*z_Flex(+1)*eta_Flex(+1);
// 20
nu_Flex=(1-theta)*betta*Lambda_Flex*(Rk_Flex(+1)-R_Flex)+betta*Lambda_Flex*theta*x_Flex(+1)*nu_Flex(+1);
// 21
phi_Flex=eta_Flex/(lambda-nu_Flex);

% Capital Producing
// 22
Q_Flex  =   1+eta_i/2*(((I_Flex-delta_Flex*exp(zeta)*K_Flex(-1))+Iss)/(In_Flex(-1)+Iss)-1)^2+eta_i*(((I_Flex-delta_Flex*exp(zeta)*K_Flex(-1))+Iss)/(In_Flex(-1)+Iss)-1)*((I_Flex-delta_Flex*exp(zeta)*K_Flex(-1))+Iss)/(In_Flex(-1)+Iss)-betta*Lambda_Flex*eta_i*((In_Flex(+1)+Iss)/((I_Flex-delta_Flex*exp(zeta)*K_Flex(-1))+Iss)-1)*((In_Flex(+1)+Iss)/((I_Flex-delta_Flex*exp(zeta)*K_Flex(-1))+Iss))^2;
// 23
Rk_Flex=(((epsilon-1)/epsilon)*alfa*(Y_Flex/K_Flex(-1))+exp(zeta)*(Q_Flex-delta_Flex))/Q_Flex(-1);
// 24
In_Flex=I_Flex-delta_Flex*exp(zeta)*K_Flex(-1);

% LoM of Capital
// 25
K_Flex=exp(zeta)*K_Flex(-1)+In_Flex;
// 26
Premium_Flex=Rk_Flex(1)-R_Flex;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sticky Prices and Sticky Wages Economy %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Households
// 27
varrho=(C-h*C(-1))^(-1)-betta*h*(C(+1)-h*C)^(-1);
// 28
Lambda=varrho(+1)/varrho;
// 29
1=betta*(varrho(+1)/varrho)*R;

% Wage Rigidities
%x1_w = (thetaw/(thetaw-1))*x2_w;
// 30
x1_w = (thetaw/(thetaw-1))*exp(markup_w)*x2_w;
// 31
x1_w= varrho*w*pop_w^(1-thetaw)*L+betta*gamma_w*((pie^(gammaW)/Pie_W(+1))*(pop_w/pop_w(+1)))^(1-thetaw)*x1_w(+1);
// 32
x2_w = chi*(pop_w)^(-thetaw*(varphi+1))*L^(1+varphi) + betta*gamma_w*((Pie_W(+1)/(pie^(gammaW)))*(pop_w(+1)/pop_w))^(thetaw*(varphi+1))*x2_w(+1);
// 33
1 = gamma_w*(Pie_W/(pie(-1)^gammaW))^(thetaw-1) + (1-gamma_w)*pop_w^(1-thetaw);
// 34
w= (Pie_W/(pie*exp(EPi)))*w(-1);


% Economy wide resource constraint
// 35
Y=C+I+(eta_i/2)*((I-delta*exp(zeta)*K(-1)+Iss)/(In(-1)+Iss)-1)^2*(I-delta*exp(zeta)*K+Iss)+tau*psi*Q*S+G;%+G; %

% % Government Expenditures
// 36
G=g*Y;
// 37
g=(1-rho_g)*g_ss+rho_g*g(-1)+sigma_epsiG*eps_g;

% Credit Policy of Central Bank
// 38
psi=psi_ss+cpp*(Premium-Premium_SS);
// 39
Q*S=phic*N;
// 40
S=Sp+psi*S;
// 41
phic=phi/(1-psi);

% Intermediate Goods-producing Firms
// 42
K=S;
%Y=exp(A)*(exp(zeta)*K(-1))^(alfa)*L^(1-alfa);   % without U
// 43
Y=exp(A)*(U*exp(zeta)*K(-1))^(alfa)*L^(1-alfa);   % with U
// 44
w=Pm*(1-alfa)*(Y/L);							% comp labor market equilibrium
// 45
Pm=exp(zeta)*K(-1)*b*U^cap*(U/(alfa*Y));
// 46
delta = delta_c+b/(1+cap)*(U^(1+cap)); % depreciation rate

% Financial Intermediaries
// 47
N=Ne+Nn;%-eps_w;
// 48
Ne=theta*z*N(-1)*exp(EWealth);
// 49
Nn=omega*Q*S(-1);
%Q*S=phi*N;
// 50
z=(Rk-R(-1))*phi(-1)+R(-1);
// 51
x=(phi/phi(-1))*z;
// 52
eta=(1-theta)+betta*Lambda*theta*z(+1)*eta(+1);
// 53
nu=(1-theta)*betta*Lambda*(Rk(+1)-R)+betta*Lambda*theta*x(+1)*nu(+1);
// 54
phi=eta/(lambda-nu);

% Capital Producing
// 55
Q  =   1+eta_i/2*(((I-delta*exp(zeta)*K(-1))+Iss)/(In(-1)+Iss)-1)^2+eta_i*(((I-delta*exp(zeta)*K(-1))+Iss)/(In(-1)+Iss)-1)*((I-delta*exp(zeta)*K(-1))+Iss)/(In(-1)+Iss)-betta*Lambda*eta_i*((In(+1)+Iss)/((I-delta*exp(zeta)*K(-1))+Iss)-1)*((In(+1)+Iss)/((I-delta*exp(zeta)*K(-1))+Iss))^2;
// 56
Rk=(Pm*alfa*(Y/K(-1))+exp(zeta)*(Q-delta))/Q(-1);
// 57
In=I-delta*exp(zeta)*K(-1);

% LoM of Capital
// 58
K=exp(zeta)*K(-1)+In;

%Retailer
%x1=(epsilon/(epsilon-1))*x2;
// 59
x1=(epsilon/(epsilon-1))*exp(markup_p)*x2;
// 60
x1=Y*(pop)^(1-epsilon)+gamma_p*betta*Lambda*(pop/pop(+1))^(1-epsilon)*(pie*exp(EPi))^(gammaP)*(pie(+1)*exp(EPi(+1)))^(epsilon-1)*x1(+1);
// 61
x2=Y*Pm*pop^(-epsilon)+gamma_p*betta*Lambda*(pie(+1)*exp(EPi(+1))*(pop(+1)/pop))^(epsilon)*x2(+1);
// 62
%1 = gamma_p*(pie*exp(EPi))^(gammaP*(1-epsilon))*(pie*exp(EPi))^(epsilon-1) + (1-gamma_p)*pop^(1-epsilon);
1 = gamma_p*(exp(EPi)*pie(-1)^gammaP/pie)^(1-epsilon) + (1-gamma_p)*(pop)^(1-epsilon);
// 63
MCInv=(1/Pm);

% Composite Inflation
// 64
pieComp=pie*exp(EPi)*(1-gamma_comp)+Pie_W*(gamma_comp);
%y_gap2=(epsilon/(epsilon-1))/MCInv;
y_gap2=Y-Y_Flex;
y_gap=Y/Y_Flex;

% Interest rate rule
i=R*pie(+1)*exp(EPi(+1));
%i=((1/betta)*pie^(kappa_pie)*(y_gap)^(kappa_y)*Pie_W^(kappa_pieW)*(Premium/steady_state(Premium))^(-kappa_prem))^(1-rho)*i(-1)^(rho)*exp(-eps_m);%*exp(-eps_m)
% Taylor rule with target of asset growth
i=((1/betta)*pie^(kappa_pie)*(y_gap)^(kappa_y)*Pie_W^(kappa_pieW)*((Rk-R(-1))/(steady_state(Rk)-steady_state(R)))^(-kappa_prem)*(x/x_Flex)^(kappa_x))^(1-rho)*i(-1)^(rho)*exp(Eint);%*exp(-eps_m)
%i=((1/betta)*(Pie_W)^(kappa_pie)*(MCInv/(epsilon/(epsilon-1)))^(-kappa_y))^(1-rho)*i(-1)^(rho)*exp(eps_m);
%R=((1/betta)*pie^(kappa_pie)*(MCInv/(epsilon/(epsilon-1)))^(-kappa_y))^(1-rho)*i(-1)^(rho)*exp(eps_m);


% Technology process
A=rho_a*A(-1)-sigma_epsiA*eps_a;
% Capital Quality AR(1)
zeta=rho_ksi*zeta(-1)-sigma_epsiK*ksi;
% ARMA(1,1) process for wage mark-up
markup_w=rho_mark*markup_w(-1)-sigma_epsiMW*eps_mark_w;%-mu_mark*eps_mark(-1);
markup_p=rho_mark_p*markup_p(-1)-sigma_epsiMP*eps_mark_p;%-mu_mark*eps_mark(-1);
% Wealth Shock AR(1)
EWealth = rho_w * EWealth(-1);% - sigma_epsiW*eps_w;
Eint = rho_Eint*Eint(-1)-sigma_epsiInt*eps_m;
EPi = rho_pie*EPi(-1);%+eps_pie;%*sigma_epsiInt*eps_m;
%EPr= rho_prem *EPr(-1)+sigma_epsiPrem*eps_prem;
% Premium
%Premium=(Rk(1)-R)*exp(EPr);
%exp(Eint) = ;
% Premium
Premium=Rk(1)-R;
%MRS=chi*L^(varphi)/varrho;
MPL=(1-alfa)*(Y/L);
MPK = alfa * (Y/K);

% Var for % deviation from St.St.
Y_log =log(Y/steady_state(Y))*100;
C_log =log(C/steady_state(C))*100;
L_log =log(L/steady_state(L))*100;
w_log =log(w/steady_state(w))*100;
Q_log =log(Q/steady_state(Q))*100;
K_log =log(K/steady_state(K))*100;
I_log =log(I/steady_state(I))*100;
N_log =log(N/steady_state(N))*100;
phic_log=log(phic/steady_state(phic))*100;
MPL_log=log(MPL/steady_state(MPL))*100;
MPK_log=log(MPK/steady_state(MPK))*100;
MC_log=log(Pm/steady_state(Pm))*100;
x_obs = log(x/steady_state(x));
Wf_obs = -log(Wf/steady_state(Wf))*100;
%Wf_obs = log(Wf);
% Welfare
Wf=log(C-h*C(-1))-chi/(1+varphi)*L^(1+varphi)+betta*Wf(+1);

EUL = chi/(1+varphi)*L^(1+varphi) +betta*EUL(+1);

C_lag=C(-1);
Incentive = (nu*(Q*S)+eta*N) - lambda*(Q*S);
Rk_nom = Rk*pie*exp(EPi);

end;
shocks;
%if sigma_epsiA != 0
    var eps_a; stderr 1;
%end
%if sigma_epsiK != 0
%    var ksi; stderr 1;
%end
%if sigma_epsiInt != 0
%    var eps_m; stderr 1;
%end
%if sigma_epsiW != 0
 %  var eps_w; stderr 1;
%end
%if sigma_epsiG != 0
 %   var eps_g; stderr 1;
%end
%if sigma_epsiMW != 0
 %  var eps_mark_w; stderr 1;
%end
%if sigma_epsiMP != 0
 %   var eps_mark_p; stderr 1;
%end
%var eps_pie; stderr 0.01;
end;

% optim_weights;
% Pie_W 1;
% pie 1;
% y_gap 1;
% end;

% osr_params kappa_pie kappa_pieW kappa_y;

%steady;
%resid;
%check;

%osr(order=1,irf=50,noprint,nodisplay) Y_log C_log L_log w_log K_log phic_log N_log pie Pie_W R Premium Q_log MPK_log MPL_log Wf;%Rk
if irf_plot == 1
    %stoch_simul(order=1,irf=50,noprint,nodisplay) Y_log C_log L_log w_log K_log phic_log N_log Rk pie Pie_W R Premium Q_log MPK_log MPL_log y_gap Wf;
    %stoch_simul(order=1,irf=40,nodisplay,noprint) Y_log C_log L_log w_log K_log phic_log N_log i pie Pie_W R Premium Q_log MPK_log MC_log Wf;
    %stoch_simul(order=1,irf=40,nodisplay,noprint) Y_log C_log L_log w_log K_log phic_log N_log i pie Pie_W R Premium Q_log MPK_log MC_log x z;
    stoch_simul(order=1,irf=40,pruning,nodisplay,noprint) Y_log C_log L_log w_log K_log phic_log Wf_obs pie Pie_W Premium R i;%y_gap
    
end
%stoch_simul(order=2,periods=10000,nodisplay,noprint) Y_log C_log L_log w_log K_log pie Pie_W R Q_log MPL_log MPK_log Wf;

if irf_plot == 0
    %stoch_simul(order=1,irf=40,pruning,nodisplay) Y_log C_log L_log w_log K_log phic_log Wf_obs pie Pie_W Premium R i;
    stoch_simul(order=2,periods=10000,nodisplay,pruning,noprint) Wf; %
    %stoch_simul(order=1,periods=1000,nodisplay) y_gap2 pie Pie_W; %
end

%stoch_simul(order=2,periods=10000,nodisplay,pruning) Wf; %
%stoch_simul(order=1,irf=20) Y_log C_log L_log w_log K_log phic_log N_log Rk pie Pie_W R Premium Q_log MPK_log MPL_log y_gap i Pm delta Rk_nom x z phi N nu eta;
%stoch_simul(order=1,irf=0,noprint,nodisplay) Y_log C_log L_log w_log K_log phic_log N_log Rk pie Pie_W R Premium i;
%stoch_simul(order=1,irf=0,noprint,nodisplay) Y_log C_log L_log w_log K_log phic_log N_log Rk pie Pie_W i Premium R Q MPK_log MPL_log;

%planner_objective(log(C-h*C_lag)-chi/(1+varphi)*L^(1+varphi));

%ramsey_policy(planner_discount=betta,order=1,instruments=(i));
