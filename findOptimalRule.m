function welfare_final = findOptimalRule(x0)

eps = 0.05;

kappa_pie   = x0(1);
kappa_pieW  = x0(2);
kappa_y     = x0(3);
kappa_prem  = x0(4);
rho         = x0(5);
kappa_x     = x0(6);

save('params','kappa_pie','kappa_pieW','kappa_x','kappa_y','kappa_prem','rho','-append');

clear global M_ oo_
%dynare GK_Nom_CapU noclearall
try
    GK_Nom_CapU
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    S_pos=strmatch('S',M_.endo_names,'exact');
    Q_pos=strmatch('Q',M_.endo_names,'exact');
    N_pos=strmatch('N',M_.endo_names,'exact');
    nu_pos=strmatch('nu',M_.endo_names,'exact');
    eta_pos=strmatch('eta',M_.endo_names,'exact');
    Incentive_pos=strmatch('Incentive',M_.endo_names,'exact');
    welfare=-(oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos)));
    S=oo_.dr.ys(S_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(S_pos));
    Q=oo_.dr.ys(Q_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(Q_pos));
    N=oo_.dr.ys(N_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(N_pos));
    nu=oo_.dr.ys(nu_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(nu_pos));
    eta=oo_.dr.ys(eta_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(eta_pos));
    
    %     S   = oo_.endo_simul(S_pos,:);
    %     Q   = oo_.endo_simul(Q_pos,:);
    %     N   = oo_.endo_simul(N_pos,:);
    %     nu  = oo_.endo_simul(nu_pos,:);
    %     eta = oo_.endo_simul(eta_pos,:);
    %     Incentive = oo_.endo_simul(Incentive_pos,:);
    %welfare=-oo_.mean(1,1);
    %iter= iter+1;
    %fprintf('Iteration: %d  \n',iter)
    param_opt_valid = x0 + eps ;
    
        kappa_pie   = param_opt_valid(1);
        kappa_pieW  = param_opt_valid(2);
        kappa_y     = param_opt_valid(3);
        %kappa_prem  = param_opt_valid(4);% 1.2
        rho         = param_opt_valid(5);
        kappa_x     = param_opt_valid(6);
    %     save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','-append');
    %     GK_Nom_CapU
    %     W_pos=strmatch('Wf',M_.endo_names,'exact');
    %     WelfOutput=-(oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos)));
    
    Term_Wealth = nu.*(Q.*S)+eta.*N;
    RHS = lambda*(Q.*S);
    
   % check_valid1 = abs(welfare/WelfOutput);
    check_valid2 = abs(Term_Wealth - RHS);
    
    %check_valid = Term_Wealth- RHS > - 0.01;
    
    %check_valid = Incentive > - 0.01;
    if check_valid2 < 0.25  %|| check_valid1 < 
        %if sum(check_valid)/length(check_valid) < 0.001
        welfare_final = welfare;
    else
        welfare_final = 100000000;
    end
catch
    welfare_final = 100000000;
end
fprintf('Welfare: %d  \n',-welfare_final)
end