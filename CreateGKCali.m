alfa        = 0.33;
betta       = 0.99;
cap         = 7.2;
chi         = 1;
cpp         = 0;
deltai      = 0.025;
epsilon     = 4.167;
eta_i       = 1.7280;
g_ss        = 0.2;
gamma_p     = 0.779; 
gamma_w     = 0.779;
gammaP      = 0.241;
gammaW      = 0.241;
weight_pie  = 0;
gamma_comp  = 0;
h           = 0.815;
i_bar       = 1/betta;
kappa_pie   = 1.5;
kappa_y     = 0.1250;
kappa_pieW  = 0;
kappa_prem  = 0;
kappa_x     = 0;
lambda      = 0.381;
%lambda      = 0.3750;
mu_mark     = 0;
omega       = 0.002;
%omega       = 0.0020;
%Premium_SS  = 150/40000;
%Premium_SS  = 0.003693602965625;
%Premium_SS  = 0.0025;
psi_ss      = 0;
pi_ss       = 0;
piW_ss      = 0;
rho         = 0.8;
rho_a       = 0.95;    
rho_Eint    = 0;    
rho_ksi     = 0.66;
rho_g       = 0;
rho_mark_p  = 0;
rho_mark    = 0;
rho_w       = 0;
rho_prem    = 0;
tau         = 1.0000e-03;
theta       = 0.972;
thetaw      = 4.167;
varphi      = 0.276;
sigma_epsiA         = 0.01;
sigma_epsiK         = 0.05;
sigma_epsiMW        = 0;
sigma_epsiPrem      = 0;
sigma_epsiMP        = 0;
sigma_epsiW         = 0;
sigma_epsiInt       = 0.0250;
sigma_epsiPrem_obs  = 0; 
sigma_epsiG         = 0;

%  phi = 4;
% 
% % lambda = ((1-theta)*betta*Premium_SS)/(1-betta*theta*(Premium_SS*phi+i_bar))+1/phi*((1-theta))/(1-betta*theta*(Premium_SS*phi+i_bar));
%  lambda = (-(1-theta)*(1+betta*Premium_SS*phi)/(betta*theta*Premium_SS*phi^2-phi*(1-theta)));
% 
%  disp(lambda)

% aa     =   lambda*betta*theta*Premium_SS;
% bb     =   -(1-theta)*(lambda-betta*Premium_SS);
% cc     =   (1-theta);
% phi    =   (-bb-sqrt(bb^2-4*aa*cc))/(2*aa);
% 
% disp(phi)

save params.mat

[ys,check,xfs]  = CapU_steadystate;
K_ss        = ys(1);
b           = ys(55);
delta_c     = ys(56);
Premium_SS  = xfs(4);
L_ss        = xfs(5);

K_ss_init   = K_ss;
Iss         = deltai*K_ss;

save GKCali alfa b betta cap chi cpp deltai delta_c epsilon eta_i g_ss gamma_p gamma_w gammaP gammaW weight_pie ...
    gamma_comp h i_bar K_ss L_ss Iss kappa_pie kappa_y kappa_pieW kappa_prem kappa_x lambda mu_mark ...
    omega Premium_SS psi_ss pi_ss piW_ss rho rho_a rho_Eint rho_ksi rho_g rho_mark_p rho_mark rho_w rho_prem tau theta thetaw varphi ...
    sigma_epsiA sigma_epsiK sigma_epsiMW sigma_epsiPrem sigma_epsiMP sigma_epsiW sigma_epsiInt sigma_epsiPrem_obs sigma_epsiG

save params.mat

ind_A     = 1;
ind_K     = 1;
ind_MW    = 1;
ind_MP    = 1;
ind_W     = 0;
ind_Int   = 1;
ind_G     = 1;

save('params.mat','ind_A','ind_K','ind_MW','ind_MP','ind_W','ind_Int','ind_G','-append')

% clear all 
% load('H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_CapU_SR\paramsBench.mat')
% %gamma_p     = 0.75;
% %gamma_w     = 0.75;
% 
% % %kappa_pie   = 1.5;
%  kappa_pieW 	= 0;
%  %kappa_y     = 0;
%  kappa_prem  = 0;
% %rho         = 0;    %
%  kappa_x     = 0;
% % %h           = 0;
% % %chi         = 1;
% % %lambda      = 0.381;
% % %omega       = 0.002;
% % %theta       = 0.975;
% % %epsilon     = 4;
% % %thetaw      = 4;
% sigma_epsiA     = 0.01;
% sigma_epsiK     = 0.01;
%  sigma_epsiMW    = 0.01;
%  sigma_epsiMP    = 0.01;
%  sigma_epsiW     = 0.01;
%  sigma_epsiInt   = 0.01;
%  sigma_epsiG       = 0.01;
%  sigma_epsiPrem  = 0.01;
%  sigma_epsiPrem_obs = 0.01;
%   pi_ss           = 0;
%   piW_ss          = 0;
%  rho_Eint       = 0;
%  rho_mark_p      = 0.9;
%  rho_mark        = 0.9;
%  rho_w           = 0.9;
%  rho_prem        = 0.9;
% %  %gammaW          = 0;
% %  %gammaP          = 0;
%   gammaP      = 0.35;
%   gammaW      = 0.39;
% save params.mat
% 
% [ys,check,xfs]  = CapU_steadystate;
% K_ss        = ys(1);
% b           = ys(55);
% delta_c     = ys(56);
% L_ss        = xfs(5);
% 
% K_ss_init   = K_ss;
% weight_pie = 0;
% 
% ys_init     = ys;
% %betta       = 0.99;
% %h           = 0;  % only small influence on the levels of IRFs, smooths consumption IRFs
% %kappa_y     = 0.5; % larger weight on output gap does not affect the sticky wages only case
% %global alfa b betta cap chi cpp deltai epsilon eta_i g_ss gamma_p gamma_w gammaP gamma_comp h i_bar K_ss Iss kappa_pie kappa_y lambda mu_mark omega Premium_SS psi_ss rho rho_a rho_ksi rho_g rho_mark rho_w tau theta thetaw varphi;
% 
% 
% %load optimal_coeff_p
% 
% save  params.mat alfa b betta cap chi cpp deltai delta_c epsilon eta_i g_ss gamma_p gamma_w gammaP gammaW weight_pie ...
%     gamma_comp h i_bar K_ss L_ss Iss kappa_pie kappa_y kappa_pieW kappa_prem kappa_x lambda mu_mark ...
%     omega Premium_SS psi_ss pi_ss piW_ss rho rho_a rho_Eint rho_ksi rho_g rho_mark_p rho_mark rho_w rho_prem tau theta thetaw varphi ...
%     sigma_epsiA sigma_epsiK sigma_epsiMW sigma_epsiPrem sigma_epsiMP sigma_epsiW sigma_epsiInt sigma_epsiPrem_obs sigma_epsiG
% 

dynare GK_Nom_CapU