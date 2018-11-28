%function diff=solveSteadyState(XX,alfa,betta,chi,cpp,delta,epsilon,eta_i,g_ss,h,i_bar,Iss,kappa_pie,kappa_y,lambda,mu_mark,omega,psi_ss,rho,rho_a,rho_ksi,rho_g,rho_mark,rho_w,tau,theta,thetaw,varphi)
function diff=solveSteadyStateCapU(XX,alfa,b,betta,cap,chi,cpp,deltai,epsilon,eta_i,g_ss,gamma_p,gamma_w,gammaP,gamma_comp,h,i_bar,K_ss,Iss,kappa_pie,kappa_y,kappa_pieW,kappa_prem,lambda,mu_mark,omega,Premium_SS,psi_ss,rho,rho_a,rho_ksi,rho_g,rho_mark,rho_w,tau,theta,thetaw,varphi)
%global alfa b betta cap chi cpp deltai epsilon eta_i g_ss gamma_p gamma_w gammaP gamma_comp h i_bar K_ss Iss kappa_pie kappa_y lambda mu_mark omega psi_ss rho rho_a rho_ksi rho_g rho_mark rho_w tau theta thetaw varphi;
K=XX(1);
L=XX(2);
Y=K^alfa*L^(1-alfa);
C=Y-deltai*K-g_ss*Y;
varrho=(1-betta*h)/((1-h)*C);
Rk=(epsilon-1)/(epsilon)*alfa*Y/K+1-deltai;
Prem=Rk-1/betta;
aa     =   lambda*betta*theta*Prem;
bb     =   -(1-theta)*(lambda-betta*Prem);
cc     =   (1-theta);
phi    =   (-bb-sqrt(bb^2-4*aa*cc))/(2*aa);
N=omega*K/(1-theta*(Prem*phi+1/betta));
%diff(1)=chi*L^varphi-(1-betta*h)/((1-h)*C)*(epsilon-1)/(epsilon)*(1-alfa)*Y/L;
diff(1)= (1-alfa)*Y/L-(thetaw/(thetaw-1))*(L^varphi/varrho)*(epsilon/(epsilon-1));
diff(2)=phi*N-K;

end