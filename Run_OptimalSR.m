clear all
close all
clear global
%oldfolder = cd('D:\users\tkoenig\GK\Estimation_All');
%oldfolder = cd('C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\Estimation_All');
%load paramsEstimate.mat
%load('C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\Calibration\VillaCali.mat')
%load('H:\ownCloud\PhD\GK_Nom_Wages\Code\Calibration\VillaCali.mat')
load('H:\ownCloud\PhD\GK_Nom_Wages\Code\Calibration\GKCali.mat')
%load paramsBench.mat
%load calibrated_params_1st
%cd(oldfolder)
%load optimal_coeff_w
%load optimal_coeff_pw_Rho_Pi
%load optimal_coeff_pw_Rho
%load optimal_coeff_pw_NoRho
load optimal_coeff_pw
%load optimal_coeff_p

%gamma_p     = 0;
%gamma_w     = 0;

% kappa_pie   = 1.0587;
% kappa_pieW 	= 23.0123;
% kappa_y     = 48.3304;
% kappa_prem  = 0;
% rho         = 0.9968;    %
% kappa_x     = 23.2302;
% % %h           = 0;
% % %chi         = 1;
% lambda      = 0.381;
% %omega       = 0.002;
% %theta       = 0.975;
% %Premium_SS  = 0.0025;
% % % %epsilon     = 4;
% % % %thetaw      = 4;
sigma_epsiA     = 0;
sigma_epsiK     = 0.01;
%  sigma_epsiMW    = 0.01;
%  sigma_epsiMP    = 0.01;
%  sigma_epsiW     = 0.01;
sigma_epsiInt   = 0;
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
% % %  %gammaW          = 0;
% % %  %gammaP          = 0;
%   gammaP      = 0;
% gammaW      = 0;
save params.mat

[ys,check,xfs]  = CapU_steadystate;
K_ss        = ys(1);
b           = ys(55);
delta_c     = ys(56);
L_ss        = xfs(5);

K_ss_init   = K_ss;
weight_pie = 0;

ys_init     = ys;
%betta       = 0.99;
%h           = 0;  % only small influence on the levels of IRFs, smooths consumption IRFs
%kappa_y     = 0.5; % larger weight on output gap does not affect the sticky wages only case
%global alfa b betta cap chi cpp deltai epsilon eta_i g_ss gamma_p gamma_w gammaP gamma_comp h i_bar K_ss Iss kappa_pie kappa_y lambda mu_mark omega Premium_SS psi_ss rho rho_a rho_ksi rho_g rho_mark rho_w tau theta thetaw varphi;


%load optimal_coeff_p

save  params.mat alfa b betta cap chi cpp deltai delta_c epsilon eta_i g_ss gamma_p gamma_w gammaP gammaW weight_pie ...
    gamma_comp h i_bar K_ss L_ss Iss kappa_pie kappa_y kappa_pieW kappa_prem kappa_x lambda mu_mark ...
    omega Premium_SS psi_ss pi_ss piW_ss rho rho_a rho_Eint rho_ksi rho_g rho_mark_p rho_mark rho_w rho_prem tau theta thetaw varphi ...
    sigma_epsiA sigma_epsiK sigma_epsiMW sigma_epsiPrem sigma_epsiMP sigma_epsiW sigma_epsiInt sigma_epsiPrem_obs sigma_epsiG


%save('params.mat','ind_A','ind_K','ind_MW','ind_MP','ind_W','ind_Int','ind_G','-append')
gridsearch  = 1;
irf_plot    = 0;
WFcost      = 0;
Robust      = 0;
comp_static = 0;
run_corr    = 0;
pol_front   = 0;
determinacy = 0;

%% Run Search for Optimal Coefficients

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('RUN GRID SEARCH!')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if gridsearch==1
    
    irf_plot = 0;
    
    
    
    for iii=1:1
        
        
        kappa_pie   = 1.5+rand(1);
        %kappa_pie   = 0;
        kappa_pieW  = 4.5+rand(1);
        kappa_y     = 0;
        %kappa_y     = randn(1) ;
        kappa_x     = 0;
        %kappa_prem  = param_opt(4);
        %kappa_prem  = 0;
        rho         = 0;
        save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','rho','kappa_prem','irf_plot','-append');
        dynare GK_Nom_CapU noclearall
        iter = 0;
        
        options = optimset('TolFun',1e-5,'Algorithm','interior-point','Display','Off');
        x0 = [kappa_pie kappa_pieW kappa_y kappa_prem rho kappa_x];
        
        kappa_pie_lb    = 1;
        kappa_pie_ub    = 50;
        %kappa_pie_lb    = 0;
        %kappa_pie_ub    = 0;
        %kappa_pieW_lb   = 0;
        %kappa_pieW_ub   = 0;
        kappa_pieW_lb   = 0;
        kappa_pieW_ub   = 50;
        kappa_y_lb      = 0;
        kappa_y_ub      = 50;
        %kappa_y_lb      = 0;
        %kappa_y_ub      = 0;
        kappa_prem_lb   = 0;
        kappa_prem_ub   = 0;
        rho_lb          = 0;
        rho_ub          = 1;
        kappa_x_lb      = 0;
        kappa_x_ub      = 50;
        %kappa_x_lb      = 0;
        %kappa_x_ub      = 0;
        
        lb = [kappa_pie_lb kappa_pieW_lb kappa_y_lb kappa_prem_lb rho_lb kappa_x_lb];
        ub = [kappa_pie_ub kappa_pieW_ub kappa_y_ub kappa_prem_ub rho_ub kappa_x_ub];
        
%         [param_opt,welfare_opt]= fmincon(@findOptimalRule,x0,[],[],[],[],lb,ub,[],options);
         
%         param_opt_vec(iii,:) = param_opt;
%         welfare_opt_vec(iii) = -welfare_opt;
        
        rng default % For reproducibility
        ms = MultiStart('UseParallel',true);
        gs = GlobalSearch(ms);
        problem = createOptimProblem('fmincon','x0',x0,'objective',@findOptimalRule,'lb',lb,'ub',ub);
        
        [param_opt,welfare_opt] = run(gs,problem);
        param_opt_vec(iii,:) = param_opt;
        welfare_opt_vec(iii) = -welfare_opt;
    end
    
    
    fprintf('Optimal Weight on Inflation: %d / \n',param_opt_vec(1,1))
    fprintf('Optimal Weight on Wage Inflation: %d / \n',param_opt_vec(1,2))
    fprintf('Optimal Weight on Output: %d / \n',param_opt_vec(1,3))
    fprintf('Optimal Weight on Premium: %d / \n',param_opt_vec(1,4))
    fprintf('Optimal Indexation: %d / \n',param_opt_vec(1,5))
    fprintf('Optimal Weight on Asset Growth: %d / \n',param_opt_vec(1,6))
    fprintf('Optimal Welfare: %d / \n',welfare_opt_vec(1))
    
%     fprintf('Optimal Weight on Inflation: %d / \n',param_opt_vec(2,1))
%     fprintf('Optimal Weight on Wage Inflation: %d / \n',param_opt_vec(2,2))
%     fprintf('Optimal Weight on Output: %d / \n',param_opt_vec(2,3))
%     fprintf('Optimal Weight on Premium: %d / \n',param_opt_vec(2,4))
%     fprintf('Optimal Indexation: %d / \n',param_opt_vec(2,5))
%     fprintf('Optimal Weight on Asset Growth: %d / \n',param_opt_vec(2,6))
%     fprintf('Optimal Welfare: %d / \n',welfare_opt_vec(2))
%     
%     fprintf('Optimal Weight on Inflation: %d / \n',param_opt_vec(3,1))
%     fprintf('Optimal Weight on Wage Inflation: %d / \n',param_opt_vec(3,2))
%     fprintf('Optimal Weight on Output: %d / \n',param_opt_vec(3,3))
%     fprintf('Optimal Weight on Premium: %d / \n',param_opt_vec(3,4))
%     fprintf('Optimal Indexation: %d / \n',param_opt_vec(3,5))
%     fprintf('Optimal Weight on Asset Growth: %d / \n',param_opt_vec(3,6))
%     fprintf('Optimal Welfare: %d / \n',welfare_opt_vec(3))
    
    %     fprintf('Optimal Weight on Inflation: %d / \n',param_opt_vec(4,1))
    %     fprintf('Optimal Weight on Wage Inflation: %d / \n',param_opt_vec(4,2))
    %     fprintf('Optimal Weight on Output: %d / \n',param_opt_vec(4,3))
    %     fprintf('Optimal Weight on Premium: %d / \n',param_opt_vec(4,4))
    %     fprintf('Optimal Indexation: %d / \n',param_opt_vec(4,5))
    %     fprintf('Optimal Weight on Asset Growth: %d / \n',param_opt_vec(4,6))
    %     fprintf('Optimal Welfare: %d / \n',welfare_opt_vec(4))
    %
    %     fprintf('Optimal Weight on Inflation: %d / \n',param_opt_vec(5,1))
    %     fprintf('Optimal Weight on Wage Inflation: %d / \n',param_opt_vec(5,2))
    %     fprintf('Optimal Weight on Output: %d / \n',param_opt_vec(5,3))
    %     fprintf('Optimal Weight on Premium: %d / \n',param_opt_vec(5,4))
    %     fprintf('Optimal Indexation: %d / \n',param_opt_vec(5,5))
    %     fprintf('Optimal Weight on Asset Growth: %d / \n',param_opt_vec(5,6))
    %     fprintf('Optimal Welfare: %d / \n',welfare_opt_vec(5))
    
    

    [M,I]=max(welfare_opt_vec(:));
    param_opt(:) = param_opt_vec(I,:);
    save optimal_coeff_pw param_opt M
    %save optimal_coeff_pw_Rho_Pi param_opt M
    
else
    disp('No grid search was running!')
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('RUN Determinancy Check!')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if determinacy==1
    
    kappa_pie_grid  = -1.5:0.125:5;
    kappa_pieW_grid = -1.5:0.125:5;
    kappa_y_grid    = -1.5:0.125:5;
    rho_grid        = 0:0.1:1;
    lambda_grid     = 0.1:0.05:0.6;
    
    indeterminacy_lambda_pie   = zeros(length(lambda_grid),length(kappa_pie_grid));
    indeterminacy_lambda_rho   = zeros(length(lambda_grid),length(rho_grid));
    indeterminacy_lambda_pieW  = zeros(length(lambda_grid),length(kappa_pieW_grid));
    indeterminacy_lambda_y     = zeros(length(lambda_grid),length(kappa_y_grid));

       
    save('params','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    iter = 0;
    
    % Determinacy Plots
    
%     for k=1:length(kappa_pie_grid)
%         kappa_pie=kappa_pie_grid(k);
%         for kk=1:length(kappa_pieW_grid)
%             kappa_pieW=kappa_pieW_grid(kk);
%             for kkk=1:length(rho_grid)
%                 rho=rho_grid(kkk);
%                 for kkkk=1:length(kappa_y_grid)
%                     kappa_y = kappa_y_grid(kkkk);
%                     save('params','kappa_pie','kappa_pieW','rho','kappa_y','-append');
%                     clear M_ oo_
%                     try
%                         GK_Nom_CapU
%                         indeterminacy(k,kk,kkk,kkkk) = 0;
%                     catch
%                         indeterminacy(k,kk,kkk,kkkk) = 1;
%                     end
%                     iter= iter+1;
%                     fprintf('Progress of Grid Search: %d / %d \n',iter,(length(kappa_pie_grid)*length(kappa_pieW_grid)*length(rho_grid))*length(kappa_y_grid))
%                    end
%             end
%         end
%     end
    
    
    for k=1:length(lambda_grid)
        lambda = lambda_grid(k);
        save('params','lambda','-append')
        [ys,check]  = CapU_steadystate;
        K_ss        = ys(1);
        b           = ys(55);
        delta_c     = ys(56);
        save('params','K_ss','b','delta_c','lambda','irf_plot','-append');

        for kk = 1:length(kappa_pie_grid)
            kappa_pie   = kappa_pie_grid(kk);
            kappa_pieW  = 0;
            kappa_y     = 0;
            kappa_prem  = 0;
            rho         = 0;
            save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','irf_plot','-append');
            try
                GK_Nom_CapU
                temp = Wf;
                indeterminacy_lambda_pie(k,kk) = 0;
            catch
                disp('Blanchard Kahn condition not fulfilled!')
                indeterminacy_lambda_pie(k,kk) = 1;
            end
            iter= iter+1;
            fprintf('Progress of Grid Search: %d / %d \n',iter,(length(kappa_pie_grid)*length(lambda_grid)))
            %clear global
        end
    end
    
    
    iter = 0;
    for k=1:length(lambda_grid)
        lambda = lambda_grid(k);
        save('params','lambda','-append')
        [ys,check]  = CapU_steadystate;
        K_ss        = ys(1);
        b           = ys(55);
        delta_c     = ys(56);
        save('params','K_ss','b','delta_c','lambda','irf_plot','-append');
        for kk = 1:length(kappa_y_grid)
            kappa_pie   = 0;
            kappa_pieW  = 0;
            kappa_y     = kappa_y_grid(kk);
            kappa_prem  = 0;
            rho         = 0;
            save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','-append');
            try
                GK_Nom_CapU
                temp = Wf;
                indeterminacy_lambda_y(k,kk) = 0;
            catch
                disp('Blanchard Kahn condition not fulfilled!')
                indeterminacy_lambda_y(k,kk) = 1;
            end
            iter= iter+1;
            fprintf('Progress of Grid Search: %d / %d \n',iter,(length(lambda_grid)*length(kappa_y_grid)))
        end
    end
    
    
    
    iter = 0;
    for k=1:length(lambda_grid)
        lambda = lambda_grid(k);
        save('params','lambda','-append')
        [ys,check]  = CapU_steadystate;
        K_ss        = ys(1);
        b           = ys(55);
        delta_c     = ys(56);
        save('params','K_ss','b','delta_c','lambda','irf_plot','-append');
        for kk = 1:length(kappa_pieW_grid)
            kappa_pie   = 0;
            kappa_pieW  = kappa_pieW_grid(kk);
            kappa_y     = 0;
            kappa_prem  = 0;
            rho         = 0;
            save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','-append');
            try
                GK_Nom_CapU
                temp = Wf;
               indeterminacy_lambda_pieW(k,kk) = 0;
            catch
                disp('Blanchard Kahn condition not fulfilled!')
                indeterminacy_lambda_pieW(k,kk) = 1;
            end
           iter= iter+1;
           fprintf('Progress of Grid Search: %d / %d \n',iter,(length(lambda_grid)*length(kappa_pieW_grid)))
         end
     end
             
    d1=figure('Name','Determinancy Lambda Pie','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,kappa_pie_grid');
    surf(X,Y,indeterminacy_lambda_pie);
    rotate3d on
    axis tight
    title('Determinancy Inflation Coefficient','interpreter','latex', 'FontSize', 16)
    xlabel('\kappa_\pi','interpreter','latex', 'FontSize', 16)
    ylabel('\lambda','interpreter','latex', 'FontSize', 16)
    xt = get(gca, 'XTick');
    %yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    %zlabel('Welfare')
    az = 0;
    el = 90;
    view(az, el);
    colormap colorcube;
    saveas(d1,'Det_Lambda_Pie','fig')
    saveas(d1,'Det_Lambda_Pie','eps')
     
    d2=figure('Name','Determinancy Lambda Y','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,kappa_y_grid');
    surf(X,Y,indeterminacy_lambda_y);
    rotate3d on
    axis tight
    title('Determinancy Output Gap Coefficient','interpreter','latex', 'FontSize', 16)
    xlabel('\kappa_\pi','interpreter','latex', 'FontSize', 16)
    ylabel('\lambda','interpreter','latex', 'FontSize', 16)
    xt = get(gca, 'XTick');
    %yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    %zlabel('Welfare')
    az = 0;
    el = 90;
    view(az, el);
    colormap colorcube;
    saveas(d2,'Det_Lambda_Y','fig')
    saveas(d2,'Det_Lambda_Y','eps')
     
    d3=figure('Name','Determinancy Lambda Pie','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,kappa_pieW_grid');
    surf(X,Y,indeterminacy_lambda_pieW);
    rotate3d on
    axis tight
    title('Determinancy Wage Inflation Coefficient','interpreter','latex', 'FontSize', 16)
    xlabel('\kappa_\pi_{W}','interpreter','latex', 'FontSize', 16)
    ylabel('\lambda','interpreter','latex', 'FontSize', 16)
    xt = get(gca, 'XTick');
    %yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    %zlabel('Welfare')
    az = 0;
    el = 90;
    view(az, el);
    colormap colorcube;
    saveas(d3,'Det_Lambda_PieW','fig')
    saveas(d3,'Det_Lambda_PieW','eps')
    
    
     % kappa_pie_grid  = -1.5:0.125:5;
    % kappa_pieW_grid = -1.5:0.125:5;
    % kappa_y_grid    = -1.5:0.125:5;
    % rho_grid        = 0:0.1:1;
    % lambda_grid     = 0.1:0.05:0.6;
    
    % indeterminacy_pie_pieW = zeros(length(kappa_pie_grid),length(kappa_pieW_grid));
    % indeterminacy_rho_pie = zeros(length(kappa_pie_grid),length(rho_grid));
    % indeterminacy_rho_pieW = zeros(length(rho_grid),length(kappa_pieW_grid));
    % indeterminacy_y_pieW= zeros(length(kappa_pieW_grid),length(kappa_y_grid));
    % indeterminacy_y_pie= zeros(length(kappa_pie_grid),length(kappa_y_grid));
    
    
    
    
    % dynare  GK_Nom_CapU
    % iter = 0;
    % lambda = 0.381;
    % [ys,check]  = CapU_steadystate;
    % K_ss        = ys(1);
    % b           = ys(55);
    % delta_c     = ys(56);
    
    % save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','irf_plot','-append');
    
    
    % for k=1:length(kappa_pieW_grid)
    %     kappa_pieW = kappa_pieW_grid(k);
    %     for kk = 1:length(kappa_pie_grid)
    %         kappa_pie   = kappa_pie_grid(kk);
    %         kappa_y     = 0;
    %         kappa_prem  = 0;
    %         rho         = 0;
    %         save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','irf_plot','-append');
    %         try
    %             clear Wf
    %             GK_Nom_CapU
    %             temp = Wf;
    %             indeterminacy_pie_pieW(k,kk) = 0;
    %         catch
    %             disp('Blanchard Kahn condition not fulfilled!')
    %             indeterminacy_pie_pieW(k,kk) = 1;
    %         end
    %         iter= iter+1;
    %         fprintf('Progress of Grid Search: %d / %d \n',iter,(length(kappa_pie_grid)*length(kappa_pieW_grid)))
    %         %clear global
    %     end
    % end
    
    % iter = 0;
    % for k=1:length(kappa_pie_grid)
    %     kappa_pie = kappa_pie_grid(k);
    %     for kk = 1:length(rho_grid)
    %         kappa_pieW  = 0;
    %         kappa_y     = 0;
    %         kappa_prem  = 0;
    %         rho         = rho_grid(kk);
    %         save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','-append');
    %         try
    %             clear Wf
    %             GK_Nom_CapU
    %             temp = Wf;
    %             indeterminacy_rho_pie(k,kk) = 0;
    %         catch
    %             disp('Blanchard Kahn condition not fulfilled!')
    %             indeterminacy_rho_pie(k,kk) = 1;
    %         end
    %         iter= iter+1;
    %         fprintf('Progress of Grid Search: %d / %d \n',iter,(length(rho_grid)*length(kappa_pie_grid)))
    %     end
    % end
    
    
    
    % iter = 0;
    % for k=1:length(rho_grid)
    %     rho = rho_grid(k);
    %     for kk = 1:length(kappa_pieW_grid)
    %         kappa_pie   = 0;
    %         kappa_pieW  = kappa_pieW_grid(kk);
    %         kappa_y     = 0;
    %         kappa_prem  = 0;
    %         save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','-append');
    %         try
    %             clear Wf
    %             GK_Nom_CapU
    %             temp = Wf;
    %             indeterminacy_rho_pieW(k,kk) = 0;
    %         catch
    %             disp('Blanchard Kahn condition not fulfilled!')
    %             indeterminacy_rho_pieW(k,kk) = 1;
    %         end
    %         iter= iter+1;
    %         fprintf('Progress of Grid Search: %d / %d \n',iter,(length(rho_grid)*length(kappa_pieW_grid)))
    %     end
    % end
    
    
    % iter = 0;
    % for k=1:length(kappa_pieW_grid)
    %     kappa_pieW = kappa_pieW_grid(k);
    %     for kk = 1:length(kappa_y_grid)
    %         kappa_pie   = 0;
    %         kappa_y     = kappa_y_grid(kk);
    %         kappa_prem  = 0;
    %         rho         = 0;
    %         save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','-append');
    %         try
    %             clear Wf
    %             GK_Nom_CapU
    %             temp = Wf;
    %             indeterminacy_y_pieW(k,kk) = 0;
    %         catch
    %             disp('Blanchard Kahn condition not fulfilled!')
    %             indeterminacy_y_pieW(k,kk) = 1;
    %         end
    %         iter= iter+1;
    %         fprintf('Progress of Grid Search: %d / %d \n',iter,(length(kappa_y_grid)*length(kappa_pieW_grid)))
    %         %clear global
    %     end
    % end
    
    
    % iter = 0;
    % for k=1:length(kappa_pie_grid)
    %     kappa_pie = kappa_pie_grid(k);
    %     for kk = 1:length(kappa_y_grid)
    %         kappa_pieW   = 0;
    %         kappa_y     = kappa_y_grid(kk);
    %         kappa_prem  = 0;
    %         rho         = 0;
    %         save('params','kappa_pie','kappa_pieW','kappa_y','lambda','rho','kappa_prem','-append');
    %         try
    %             clear Wf
    %             GK_Nom_CapU
    %             temp = Wf;
    %             indeterminacy_y_pie(k,kk) = 0;
    %         catch
    %             disp('Blanchard Kahn condition not fulfilled!')
    %             indeterminacy_y_pie(k,kk) = 1;
    %         end
    %         iter= iter+1;
    %         fprintf('Progress of Grid Search: %d / %d \n',iter,(length(kappa_y_grid)*length(kappa_pie_grid)))
    %         %clear global
    %     end
    % end
    
    % d1=figure('Name','Determinancy Pie PieW','NumberTitle','off');
    % [X,Y] = meshgrid(kappa_pie_grid,kappa_pieW_grid');
    % surf(X,Y,indeterminacy_pie_pieW);
    % rotate3d on
    % axis tight
    % title('Determinancy Price Inflation Coefficient and Wage Inflation Coefficient','interpreter','latex', 'FontSize', 16)
    % xlabel('\kappa_\pi','interpreter','latex', 'FontSize', 16)
    % ylabel('\kappa_{\pi_w}','interpreter','latex', 'FontSize', 16)
    % xt = get(gca, 'XTick');
    % %yt = get(gca, 'YTick');
    % set(gca, 'FontSize', 16)
    % %zlabel('Welfare')
    % az = 0;
    % el = 90;
    % view(az, el);
    % colormap colorcube;
    % saveas(d1,'Det_Pie_PieW','fig')
    % saveas(d1,'Det_Pie_PieW','eps')
    
    % d2=figure('Name','Determinancy Pie Rho','NumberTitle','off');
    % [X,Y] = meshgrid(kappa_pie_grid,rho_grid);
    % surf(X,Y,indeterminacy_rho_pie');
    % rotate3d on
    % axis tight
    % title('Determinancy Price Inflation Coefficient and Interest Rate Smoothing','interpreter','latex', 'FontSize', 16)
    % xlabel('\kappa_\pi','interpreter','latex', 'FontSize', 16)
    % ylabel('\rho','interpreter','latex', 'FontSize', 16)
    % xt = get(gca, 'XTick');
    % %yt = get(gca, 'YTick');
    % set(gca, 'FontSize', 16)
    % %zlabel('Welfare')
    % az = 0;
    % el = 90;
    % view(az, el);
    % colormap colorcube;
    % saveas(d2,'Det_Pie_Rho','fig')
    % saveas(d2,'Det_Pie_Rho','eps')
    
    % d3=figure('Name','Determinancy PieW Rho','NumberTitle','off');
    % [X,Y] = meshgrid(kappa_pieW_grid,rho_grid);
    % surf(X,Y,indeterminacy_rho_pieW);
    % rotate3d on
    % axis tight
    % title('Determinancy Wage Inflation Coefficient and Interest Rate Smoothing','interpreter','latex', 'FontSize', 16)
    % xlabel('\kappa_{\pi_w}','interpreter','latex', 'FontSize', 16)
    % ylabel('\rho','interpreter','latex', 'FontSize', 16)
    % xt = get(gca, 'XTick');
    % %yt = get(gca, 'YTick');
    % set(gca, 'FontSize', 16)
    % %zlabel('Welfare')
    % az = 0;
    % el = 90;
    % view(az, el);
    % colormap colorcube;
    % saveas(d3,'Det_PieW_Rho','fig')
    % saveas(d3,'Det_PieW_Rho','eps')
    
    % d4=figure('Name','Determinancy PieW Output','NumberTitle','off');
    % [X,Y] = meshgrid(kappa_pieW_grid,kappa_y_grid);
    % surf(X,Y,indeterminacy_y_pieW');
    % rotate3d on
    % axis tight
    % title('Determinancy Wage Inflation Coefficient and Output Gap Coefficient','interpreter','latex', 'FontSize', 16)
    % xlabel('\kappa_{\pi_w}','interpreter','latex', 'FontSize', 16)
    % ylabel('\kappa_y','interpreter','latex', 'FontSize', 16)
    % xt = get(gca, 'XTick');
    % %yt = get(gca, 'YTick');
    % set(gca, 'FontSize', 16)
    % %zlabel('Welfare')
    % az = 0;
    % el = 90;
    % view(az, el);
    % colormap colorcube;
    % saveas(d4,'Det_PieW_Y','fig')
    % saveas(d4,'Det_PieW_Y','eps')
    
    % d5=figure('Name','Determinancy Pie Output','NumberTitle','off');
    % [X,Y] = meshgrid(kappa_pie_grid,kappa_y_grid);
    % surf(X,Y,indeterminacy_y_pie');
    % rotate3d on
    % axis tight
    % title('Determinancy Price Inflation Coefficient and Output Gap Coefficient','interpreter','latex', 'FontSize', 16)
    % xlabel('\kappa_{\pi}','interpreter','latex', 'FontSize', 16)
    % ylabel('\kappa_y','interpreter','latex', 'FontSize', 16)
    % xt = get(gca, 'XTick');
    % %yt = get(gca, 'YTick');
    % set(gca, 'FontSize', 16)
    % %zlabel('Welfare')
    % az = 0;
    % el = 90;
    % view(az, el);
    % colormap colorcube;
    % saveas(d5,'Det_Pie_Y','fig')
    % saveas(d5,'Det_Pie_Y','eps')
    
else
    disp('No Determinancy Calculations')
end

%% Run Impulse Responses

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Impulse Responses')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

if irf_plot==1
    irf=classIRF;
    irf.name={'Output','Consumption','Employment','Real Wage','Capital','Leverage Ratio','Welfare','Inflation','Wage Inflation','Spread','Real Rate','Policy Rate'};%,'Return of Capital'
    %irf.name={'Output','Consumption','Employment','Real Wage','Capital','Leverage Ratio','Net Worth','Policy','Inflation','Wage Inflation','Real Rate','Premium','Price of Capital','MPK','MC','Welfare'};%,'Return of Capital'
    irf.shocks=1;
    irf.periods=40;
    irf.numVar=length(irf.name);
    irf.sim=cell(irf.shocks*length(irf.name),3);%*length(varphi_s));
    irf.steadyState=cell(56,3);
    oldfolder=cd('Figures');
    cd(oldfolder);
    
    load optimal_coeff_pw
    kappa_pie   = param_opt(1);
    kappa_pieW  = param_opt(2);
    kappa_y     = param_opt(3);
    kappa_prem  = param_opt(4);% 1.2
    rho         = param_opt(5);
    kappa_x     = param_opt(6);
    
%     kappa_pie   = 10; % 1.5
%     kappa_pieW  = 0;
%     kappa_y     = 0; %0.125
%     kappa_prem  = 0;% 1.2
%     rho         = 0; % 0.8
%     kappa_x       = 0;
%     gamma_w =0;
%     gammaW =0;
    irf_plot = 1;
    save('params','gamma_w','gammaW','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,1)=struct2cell(oo_.irfs);
    
    
    %oldfolder=cd('D:\GK\GK_CapU_Flex');
    %clear params
    load params
    kappa_pie   = 1.5; % 1.5
    kappa_pieW  = 0;
    kappa_y     = 0.125; %0.125
    kappa_prem  = 0;% 1.2
    rho         = 0.8; % 0.8
    kappa_x       = 0;
    irf_plot = 1;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,2)=struct2cell(oo_.irfs);
    
    
    %cd('D:\GK\GK_CapU_Flex')
    %clear params
    load params
    kappa_pie   = 1.71;
    kappa_pieW  = 0;
    kappa_y     = 0;
    kappa_prem  = 0;% 1.2
    rho         = 0;
    kappa_x     = 0;
    irf_plot = 1;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,3)=struct2cell(oo_.irfs);
    
    
    load params
    kappa_pie   = 0;
    kappa_pieW  = 1.5; % 1.5
    kappa_y     = 0;
    kappa_prem  = 0;% 1.2
    rho         = 0;
    kappa_x       = 0;
    irf_plot = 1;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,4)=struct2cell(oo_.irfs);
    
    
    cd('Figures\Sticky Wages and Sticky Prices');
    
    % Define properties of plot:
    spaceH=0.03;
    spaceV=0.05;
    marTop=0.03;
    marBot=0.1;
    marginL=0.1;
    margin=0.1;
    padding=0.0;
    %margin=0.03;
    %marginL=0.075;
    
    h1=figure('Name','Capital Quality Shock Bench','NumberTitle','off');
    for j=1:length(irf.name)
        subaxis(4,3,j,'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
        hold on
        %for jj=1:(length(gamma_ws)*length(gamma_s))
        for jj=1:4
            if jj==3
                h3=plot(1:irf.periods,irf.sim{j,jj},'b:','LineWidth',2);
            elseif jj==1
                h1=plot(1:irf.periods,irf.sim{j,jj},'r','LineWidth',2);
            elseif jj==4
                h4=plot(1:irf.periods,irf.sim{j,jj},'k-.','LineWidth',2);
            else
                h2=plot(1:irf.periods,irf.sim{j,jj},'m--','LineWidth',2);
            end
            line([0 irf.periods],[0 0],'Color','k')
        end
        axis tight
        hold off
        if j >  length(irf.name)-3
            xlabel('quarters','interpreter','latex', 'FontSize', 16)
        end
        set(gca,'fontsize',6)
        if j<8
            ylabel('$ \% $dev.','interpreter','latex', 'FontSize', 12)
        else
            ylabel('pp. dev.','interpreter','latex', 'FontSize', 12)
        end
        title(irf.name(j),'interpreter','latex', 'FontSize', 16)
        if j == length(irf.name)-1
            l=legend([h1 h2 h3 h4],'Optimal SR','Benchmark','Inflation Targeting','Wage Targeting','Orientation','horizontal');
            set(l,'interpreter','latex','FontSize',18)
        end
    end
    saveas(h1,'GK_ksi_irf','fig')
    saveas(h1,'GK_ksi_irf','jpg')
    saveas(h1,'GK_ksi_irf','png')
    saveas(h1,'GK_ksi_irf','eps')
    

    % h2=figure('Name','Capital Quality Shock Bench','NumberTitle','off');
    % for j=length(irf.name)+1:2*length(irf.name)
    % 	subaxis(5,3,j-length(irf.name),'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
    % 	hold on
    % 	for jj=1:(length(gamma_ws)*length(gamma_s))
    % 		if jj==3
    %  			plot(1:irf.periods,irf.sim{j,jj},'b:','LineWidth',1.5)
    %  		elseif jj==4
    %  			plot(1:irf.periods,irf.sim{j,jj},'r','LineWidth',1.5)
    %  		elseif jj==1
    %  			plot(1:irf.periods,irf.sim{j,jj},'b--','LineWidth',1.5)
    %  		else
    %  			plot(1:irf.periods,irf.sim{j,jj},'r--')
    %  		end;
    %  		line([0 irf.periods],[0 0],'Color','k')
    % 	end
    % 	axis tight
    % 	hold off
    % 	%xlabel('quarters','interpreter','latex', 'FontSize', 12)
    % 	set(gca,'fontsize',6)
    % 	if j<length(irf.name)+12
    %  		ylabel('$ \% $dev. from St.St','interpreter','latex', 'FontSize', 8)
    %  	else
    %  		ylabel('dev. from St.St','interpreter','latex', 'FontSize', 8)
    %  	end;
    %  	title(irf.name(j-length(irf.name)),'interpreter','latex', 'FontSize', 8)
    %  	%l=legend('$\gamma_w=0$','$\gamma_w=0.8$');
    %  	%set(l,'Interpreter','latex','FontSize',10);
    %  	%l=legend('$\gamma_w1=0$','$\gamma_w1=0.8$','$\gamma_w2=0$','$\gamma_w2=0.8$','$\gamma_w3=0$','$\gamma_w3=0.8$');
    %  	%set(l,'interpreter','latex','FontSize',10);
    %  end
    % saveas(h2,'GK_gov_ksi_Bench','fig')
    % saveas(h2,'GK_gov_ksi_Bench','png')
    
    % h3=figure('Name','Monetary Policy Shock Bench','NumberTitle','off');
    % for j=2*length(irf.name)+1:3*length(irf.name)
    % 	subaxis(5,3,j-2*length(irf.name),'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
    % 	hold on
    % 	for jj=1:(length(gamma_ws)*length(gamma_s))
    % 		if jj==3
    %  			plot(1:irf.periods,irf.sim{j,jj},'b:','LineWidth',2)
    %  		elseif jj==4
    %  			plot(1:irf.periods,irf.sim{j,jj},'r','LineWidth',2)
    %  		elseif jj==1
    %  			plot(1:irf.periods,irf.sim{j,jj},'b--','LineWidth',2)
    %  		else
    %  			plot(1:irf.periods,irf.sim{j,jj},'r--')
    %  		end;
    %  		line([0 irf.periods],[0 0],'Color','k')
    % 	end
    % 	axis tight
    % 	hold off
    % 	%xlabel('quarters','interpreter','latex', 'FontSize', 12)
    % 	set(gca,'fontsize',6)
    % 	if j<2*length(irf.name)+12
    %  		ylabel('$ \% $dev. from St.St','interpreter','latex', 'FontSize', 8)
    %  	else
    %  		ylabel('dev. from St.St','interpreter','latex', 'FontSize', 8)
    %  	end;
    %  	title(irf.name(j-2*length(irf.name)),'interpreter','latex', 'FontSize', 8)
    %  	%l=legend('$\gamma_w=0$','$\gamma_w=0.8$');
    %  	%set(l,'Interpreter','latex','FontSize',10);
    %  	%l=legend('$\gamma_w1=0$','$\gamma_w1=0.8$','$\gamma_w2=0$','$\gamma_w2=0.8$','$\gamma_w3=0$','$\gamma_w3=0.8$');
    %  	%set(l,'interpreter','latex','FontSize',10);
    %  end
    % saveas(h3,'GK_gov_eps_m_Bench','fig')
    % saveas(h3,'GK_gov_eps_m_Bench','jpg')
    
    % h4=figure('Name','Wealth Shock Bench','NumberTitle','off');
    % for j=3*length(irf.name)+1:4*length(irf.name)
    % 	subaxis(5,3,j-3*length(irf.name),'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
    % 	hold on
    % 	for jj=1:(length(gamma_ws)*length(gamma_s))
    % 		if jj==3
    %  			plot(1:irf.periods,irf.sim{j,jj},'b:','LineWidth',2)
    %  		elseif jj==4
    %  			plot(1:irf.periods,irf.sim{j,jj},'r','LineWidth',2)
    %  		elseif jj==1
    %  			plot(1:irf.periods,irf.sim{j,jj},'b--','LineWidth',2)
    %  		else
    %  			plot(1:irf.periods,irf.sim{j,jj},'r--')
    %  		end;
    %  		line([0 irf.periods],[0 0],'Color','k')
    % 	end
    % 	axis tight
    % 	hold off
    % 	%xlabel('quarters','interpreter','latex', 'FontSize', 12)
    % 	if j<3*length(irf.name)+10
    %  		ylabel('$ \% $dev. from St.St','interpreter','latex', 'FontSize', 10)
    %  	else
    %  		ylabel('dev. from St.St','interpreter','latex', 'FontSize', 10)
    %  	end;
    %  	title(irf.name(j-3*length(irf.name)),'interpreter','latex', 'FontSize', 12)
    %  	%l=legend('$\gamma_w=0$','$\gamma_w=0.8$');
    %  	%set(l,'Interpreter','latex','FontSize',10);
    %  	%l=legend('$\gamma_w1=0$','$\gamma_w1=0.8$','$\gamma_w2=0$','$\gamma_w2=0.8$','$\gamma_w3=0$','$\gamma_w3=0.8$');
    %  	%set(l,'interpreter','latex','FontSize',10);
    %  end
    % saveas(h4,'GK_gov_eps_w','fig')
    % saveas(h4,'GK_gov_eps_w','jpg')
    
    % h5=figure('Name','Wage Mark-up Shock Bench','NumberTitle','off');
    % for j=4*length(irf.name)+1:5*length(irf.name)
    % 	subaxis(5,3,j-4*length(irf.name),'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
    % 	hold on
    % 	for jj=1:(length(gamma_ws)*length(gamma_s))
    % 		if jj==3
    %  			plot(1:irf.periods,irf.sim{j,jj},'b:','LineWidth',2)
    %  		elseif jj==4
    %  			plot(1:irf.periods,irf.sim{j,jj},'r','LineWidth',2)
    %  		elseif jj==1
    %  			plot(1:irf.periods,irf.sim{j,jj},'b--','LineWidth',2)
    %  		else
    %  			plot(1:irf.periods,irf.sim{j,jj},'r--')
    %  		end;
    %  		line([0 irf.periods],[0 0],'Color','k')
    % 	end
    % 	axis tight
    % 	hold off
    % 	xlabel('quarters','interpreter','latex', 'FontSize', 12)
    % 	if j<4*length(irf.name)+10
    %  		ylabel('$ \% $dev. from St.St','interpreter','latex', 'FontSize', 10)
    %  	else
    %  		ylabel('dev. from St.St','interpreter','latex', 'FontSize', 10)
    %  	end;
    %  	title(irf.name(j-4*length(irf.name)),'interpreter','latex', 'FontSize', 12)
    %  	%l=legend('$\gamma_w=0$','$\gamma_w=0.8$');
    %  	%set(l,'Interpreter','latex','FontSize',10);
    %  	%l=legend('$\gamma_w1=0$','$\gamma_w1=0.8$','$\gamma_w2=0$','$\gamma_w2=0.8$','$\gamma_w3=0$','$\gamma_w3=0.8$');
    %  	%set(l,'interpreter','latex','FontSize',10);
    %  end
    % saveas(h5,'GK_gov_eps_markup','fig')
    % saveas(h5,'GK_gov_eps_markup','jpg')
    
    cd(oldfolder)
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Nominal Rigidities %
    %%%%%%%%%%%%%%%%%%%%%%
    
    irf_noFF = irf.sim;
    
    load params
    kappa_pie   = 1.5;
    kappa_pieW  = 0;
    kappa_y     = 0;
    kappa_prem  = 0;% 1.2
    rho         = 0;
    kappa_x       = 0;
    irf_plot = 1;
    gamma_p     = 0.75;
    gamma_w     = 0.75;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','gamma_p','gamma_w','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,1)=struct2cell(oo_.irfs);
    %copyfile params.mat H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    copyfile params.mat C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    %oldfolder=cd('H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    oldfolder=cd('C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    clear oo_ M_
    [ys,~]  = CapU_steadystate;
    K_ss        = ys(1);
    b           = ys(33);
    delta_c     = ys(34);
    irf_plot = 1;
    save('params','K_ss','b','delta_c','irf_plot','-append');
    dynare NK_NoFinFric noclearall
    irf_noFF(:,1) =struct2cell(oo_.irfs);
    cd(oldfolder)
    
    
    %cd('D:\GK\GK_CapU_Flex')
    %clear params
    % No nom. rigidities
    load params
    gamma_p     = 0;
    gamma_w     = 0;
    irf_plot = 1;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','gamma_p','gamma_w','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,2)=struct2cell(oo_.irfs);
    %copyfile params.mat H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    copyfile params.mat C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    %oldfolder=cd('H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    oldfolder=cd('C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    clear oo_ M_
    [ys,~]  = CapU_steadystate;
    K_ss        = ys(1);
    b           = ys(33);
    delta_c     = ys(34);
    irf_plot = 1;
    save('params','K_ss','b','delta_c','irf_plot','-append');
    dynare NK_NoFinFric noclearall
    irf_noFF(:,2) =struct2cell(oo_.irfs);
    cd(oldfolder)
    
    % Only sticky prices
    load params
    gamma_p     = 0.75;
    gamma_w     = 0;
    irf_plot = 1;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','gamma_p','gamma_w','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,3)=struct2cell(oo_.irfs);
   % copyfile params.mat H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    copyfile params.mat C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    %oldfolder=cd('H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    oldfolder=cd('C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    clear oo_ M_
    [ys,~]  = CapU_steadystate;
    K_ss        = ys(1);
    b           = ys(33);
    delta_c     = ys(34);
    irf_plot = 1;
    save('params','K_ss','b','delta_c','irf_plot','-append');
    dynare NK_NoFinFric noclearall
    irf_noFF(:,3) =struct2cell(oo_.irfs);
    cd(oldfolder)
    
    % Only sticky wages 
    load params
    gamma_p     = 0;
    gamma_w     = 0.75;
    irf_plot = 1;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','gamma_p','gamma_w','-append');
    dynare GK_Nom_CapU noclearall
    irf.sim(:,4)=struct2cell(oo_.irfs);
    %copyfile params.mat H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    copyfile params.mat C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw
    %oldfolder=cd('H:\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    oldfolder=cd('C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    clear oo_ M_
    [ys,~]  = CapU_steadystate;
    K_ss        = ys(1);
    b           = ys(33);
    delta_c     = ys(34);
    irf_plot = 1;
    save('params','K_ss','b','delta_c','irf_plot','-append');
    dynare NK_NoFinFric noclearall
    irf_noFF(:,4) =struct2cell(oo_.irfs);
    cd(oldfolder)
    
    cd('Figures\Sticky Wages and Sticky Prices');
    
    % Define properties of plot:
    spaceH=0.03;
    spaceV=0.05;
    marTop=0.03;
    marBot=0.1;
    marginL=0.1;
    margin=0.1;
    padding=0.0;
    %margin=0.03;
    %marginL=0.075;
    
    h1=figure('Name','Monetary shock with different rigidities','NumberTitle','off');
    for j=1:length(irf.name)
        subaxis(4,3,j,'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
        hold on
        %for jj=1:(length(gamma_ws)*length(gamma_s))
        for jj=1:4
            if jj==3
                h3=plot(1:irf.periods,irf.sim{j,jj},'b:','LineWidth',2);
            elseif jj==1
                h1=plot(1:irf.periods,irf.sim{j,jj},'r','LineWidth',2);
            elseif jj==4
                h4=plot(1:irf.periods,irf.sim{j,jj},'k-.','LineWidth',2);
            else
                h2=plot(1:irf.periods,irf.sim{j,jj},'m--','LineWidth',2);
            end
            line([0 irf.periods],[0 0],'Color','k')
        end
        axis tight
        hold off
        if j >  length(irf.name)-3
            xlabel('quarters','FontSize', 14)
        end
        set(gca,'fontsize',6)
        if j<8
            ylabel('$ \% $dev.','interpreter','latex', 'FontSize', 12)
        else
            ylabel('pp. dev.','interpreter','latex', 'FontSize', 12)
        end
        title(irf.name(j),'interpreter','latex', 'FontSize', 12)
        if j == length(irf.name)-1
            l=legend([h1 h2 h3 h4],'both nom. rigid.','No nom. rigid.','only prices','only wages','Orientation','horizontal');
            set(l,'interpreter','latex','FontSize',14)
        end
    end
    saveas(h1,'GK_ksi_irf_nom','fig')
    saveas(h1,'GK_ksi_irf_nom','jpg')
    saveas(h1,'GK_ksi_irf_nom','png')
    saveas(h1,'GK_ksi_irf_nom','eps')
    
    
h1=figure('Name','Monetary Shock with and without FF','NumberTitle','off');
    for j=1:length(irf.name)
        subaxis(4,3,j,'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
        hold on
        %for jj=1:(length(gamma_ws)*length(gamma_s))
        for jj=1:4
            % Sticky prices only
            if jj==3
                h3=plot(1:irf.periods,irf.sim{j,jj},'b:','LineWidth',2);
                h1=plot(1:irf.periods,irf_noFF{j,jj},'r','LineWidth',2);
                         
             % Sticky wages only 
             elseif jj==4
                 h4=plot(1:irf.periods,irf.sim{j,jj},'k-.','LineWidth',2);
                 h2=plot(1:irf.periods,irf_noFF{j,jj},'m--','LineWidth',2);
            end
            line([0 irf.periods],[0 0],'Color','k')
        end
        axis tight
        hold off
        if j >  length(irf.name)-3
            xlabel('quarters','FontSize', 14)
        end
        set(gca,'fontsize',6)
        if j<8
            ylabel('$ \% $dev.','interpreter','latex', 'FontSize', 12)
        else
            ylabel('pp. dev.','interpreter','latex', 'FontSize', 12)
        end
        title(irf.name(j),'interpreter','latex', 'FontSize', 12)
        if j == length(irf.name)-1
            l=legend([h1 h2 h3 h4],'prices, no FF','wages, no FF','prices, FF','wages, FF','Orientation','horizontal');
            set(l,'interpreter','latex','FontSize',14)
        end
    end
    saveas(h1,'GK_ksi_irf_nomFF','fig')
    saveas(h1,'GK_ksi_irf_nomFF','jpg')
    saveas(h1,'GK_ksi_irf_nomFF','png')
    saveas(h1,'GK_ksi_irf_nomFF','eps')
    
    
    cd(oldfolder)
else
    disp('No IRFs are printed')
end

%% Run Welfare Cost Computation under Different Policies

if WFcost==1
    
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Calculate Welfare Costs Under Different Policies!')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    WelfOutput =zeros(1,13);
    
    irf_plot = 0;
    
    
    kappa_pie   = 1.5;
    kappa_y     = 0.125;
    kappa_pieW  = 0;
    kappa_prem  = 0;
    rho         = 0.8;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,1)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,1)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    
    kappa_pie   = 1.5;
    kappa_y     = 0;
    kappa_pieW  = 0;
    kappa_prem  = 0;
    rho         = 0;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,2)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,2)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    kappa_pie   = 1.5;
    kappa_y     = 0;
    kappa_pieW  = 0;
    kappa_prem  = 0;
    rho         = 0.8;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,3)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,3)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    kappa_pie   = 0;
    kappa_y     = 0;
    kappa_pieW  = 1.5;
    kappa_prem  = 0;
    rho         = 0;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,4)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,4)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    
    kappa_pie   = 0;
    kappa_y     = 0;
    kappa_pieW  = 1.5;
    kappa_prem  = 0;
    rho         = 0.8;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,5)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,5)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    kappa_pie   = 0;
    kappa_y     = 0.125;
    kappa_pieW  = 1.5;
    kappa_prem  = 0;
    rho         = 0.8;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,6)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,6)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    
    kappa_pie   = 1.5;
    kappa_y     = 0;
    kappa_pieW  = 1.5;
    kappa_prem  = 0;
    rho         = 0.8;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,7)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,7)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    
    kappa_pie   = 1;
    kappa_y     = 5;
    kappa_pieW  = 0;
    kappa_prem  = 0;
    rho         = 0;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,8)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,8)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    kappa_pie   = 0;
    kappa_y     = 0;
    kappa_pieW  = 5;
    kappa_prem  = 0;
    rho         =  0;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,9)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,9)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    

    
    
    
    kappa_pie   = param_opt(1);
    %kappa_pie   = 0;
    kappa_pieW  = param_opt(2);
    kappa_y     = param_opt(3);
    kappa_prem  = param_opt(4);% 1.2
    rho         = param_opt(5);
    kappa_x       = param_opt(6);
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,10)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,10)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    kappa_pie   = 1.5;
    kappa_y     = 0;
    kappa_pieW  = 0;
    kappa_prem  = 0;
    rho         = 0;
    kappa_x     = 0.5;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,2)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,11)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    kappa_pie   = 0;
    kappa_y     = 0;
    kappa_pieW  = 1.5;
    kappa_prem  = 0;
    rho         = 0;
    kappa_x     = 0.5;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,2)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,12)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    kappa_pie   = 5;
    kappa_y     = 0;
    kappa_pieW  = 0;
    kappa_prem  = 0;
    rho         =  0;
    kappa_x       = 0;
    save('params','K_ss','b','delta_c','kappa_x','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    %WelfOutput(1,9)=oo_.mean(1,1);
    W_pos=strmatch('Wf',M_.endo_names,'exact');
    WelfOutput(1,13)=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    
    
    %     oldfolder=cd('D:\GK\GK_CapU_Ramsey')
    %     save  params.mat alfa b betta cap chi cpp deltai epsilon eta_i g_ss gamma_p gamma_w gammaP ...
    %         gamma_comp h i_bar K_ss Iss kappa_pie kappa_y kappa_pieW kappa_prem lambda mu_mark omega Premium_SS psi_ss rho rho_a rho_ksi rho_g rho_mark rho_w tau theta thetaw varphi
    %
    %     dynare GK_Nom_CapU noclearall
    %     WfRamsey = oo_.planner_objective_value(2);
    %     cd(oldfolder)
    %     kappa_pie   = 0;
    %     kappa_y     = 0;
    %     kappa_pieW  = 0;
    %     kappa_prem  = 1.1;
    %     rho         = 0;
    %     save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','-append');
    %     dynare GK_Nom_CapU %noclearall
    %     WelfOutput(1,8)=oo_.mean(1,1);
    
    
    
    
    
    
    
    Model = {'GK Bench with Smoothing';'Price Inflation Only';'Price Inflation w/ Smoothing';'Wage Inflation Only';'Wage Inflation w/ Smoothing';'Wage Inflation, Output and Smoothing';'Wage Inflation, Price Inflation and Smoothing';'Strict Output Gap Targeting';'Strict Wage Targeting';'Optimal Simpel Rule';'Price Inflation w/ Asset Growth';'Wage Inflation w/ Asset Growth';'Strict Inflation'};
    welfare_table  = struct(    'Bench',                             [WelfOutput(:,1)],...
        'Pie',                 [WelfOutput(:,2)],...
        'PieW',                  [WelfOutput(:,3)],...
        'PieWandOutput',            [WelfOutput(:,4)],...
        'PieWandPie',   [WelfOutput(:,5)])
    
    for j=1:size(WelfOutput,2)
        WelfCost(j) =100*(1-exp((WelfOutput(1,j)-WelfOutput(1,10))*(1-betta)));
    end;
    
    WelfareCosts =WelfCost';
    T = table(WelfareCosts,'RowNames',Model)
    
    fprintf('Consumption equivalent: %d \n',WelfCost)
    
    
    %pause;
else
    disp('No welfare costs are calculated!')
    
end

%% Robustness Checks
if Robust == 1
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Calculate Robustness Checks!')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    irf_plot = 0;
    
    steps = 10;
    
    kappa_pie_grid  = 0:50/steps:50;
    kappa_pieW_grid = 0:50/steps:50;
    kappa_y_grid    = 0:50/steps:50;
    kappa_prem_grid = 0:0.05:1;
    rho_grid        = 0:1/steps:1;
    kappa_x_grid    = 0:50/steps:50;
    
    WelfCostPie     = zeros(length(kappa_pie_grid),1);
    WelfCostPieW    = zeros(length(kappa_pieW_grid),1);
    WelfCostPie2     = zeros(length(kappa_pie_grid),1);
    WelfCostPieW2   = zeros(length(kappa_pieW_grid),1);
    WelfCostY       = zeros(length(kappa_y_grid),1);
    WelfCostPrem    = zeros(length(kappa_prem_grid),1);
    WelfCostRho     = zeros(length(rho_grid),1);
    WelfCostKappx     = zeros(length(kappa_x_grid),1);
    
    %welfare_grid=zeros(length(kappa_pie_grid),length(kappa_pieW_grid),length(kappa_y_grid),length(kappa_prem_grid),length(rho_grid));
    WelfCostFig2=zeros(length(kappa_pie_grid),length(kappa_pieW_grid),length(kappa_y_grid),length(kappa_prem_grid),length(rho_grid));
    
    [ys,check]  = CapU_steadystate;
    K_ss        = ys(1);
    b           = ys(55);
    delta_c     = ys(56);
    kappa_pieW  = 1.5;
    kappa_pie   = 1.5;
    kappa_y     = 0.125;
    kappa_prem  = 0;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','irf_plot','-append');
    dynare GK_Nom_CapU noclearall
    iter = 0;
    
    for k=1:length(kappa_pie_grid)
        kappa_pie   = kappa_pie_grid(k);
        kappa_pieW  = param_opt(2);
        kappa_y     = param_opt(3);
        kappa_prem  = param_opt(4);
        rho         = param_opt(5);
        kappa_x       = param_opt(6);
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            %welfare_grid=oo_.mean(1,1);
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
            WelfCostPie(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
            WelfCostPie(k) = NaN;
        end
    end
    
        for k=1:length(kappa_pie_grid)
        kappa_pie   = kappa_pie_grid(k);
        kappa_pieW  = 0;
        kappa_y     = 0;
        kappa_prem  = 0;
        rho         = 0;
        kappa_x     = 0;
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            %welfare_grid=oo_.mean(1,1);
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
            WelfCostPie2(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
            WelfCostPie2(k) = NaN;
        end
        end
    
                for k=1:length(kappa_pie_grid)
        kappa_pie   = 0;
        kappa_pieW  = kappa_pieW_grid(k);
        kappa_y     = 0;
        kappa_prem  = 0;
        rho         = 0;
        kappa_x     = 0;
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            %welfare_grid=oo_.mean(1,1);
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
             WelfCostPieW2(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
             WelfCostPieW2(k) = NaN;
        end
    end
    
    for k=1:length(kappa_pieW_grid)
        kappa_pieW  = kappa_pieW_grid(k);
        kappa_pie   = param_opt(1);
        kappa_y     = param_opt(3);
        kappa_prem  = param_opt(4);
        rho         = param_opt(5);
        kappa_x       = param_opt(6);
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
            WelfCostPieW(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
            WelfCostPieW(k) = NaN;
        end
    end
    
    
    for k=1:length(kappa_y_grid)
        kappa_y     =kappa_y_grid(k);
        kappa_pie   = param_opt(1);
        kappa_pieW  = param_opt(2);
        kappa_prem  = param_opt(4);
        rho         = param_opt(5);
        kappa_x       = param_opt(6);
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
            WelfCostY(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
            WelfCostY(k) = NaN;
        end
    end
    
    
    
    for k=1:length(kappa_prem_grid)
        kappa_prem  = kappa_prem_grid(k);
        kappa_pie   = param_opt(1);
        kappa_pieW  = param_opt(2);
        kappa_y     = param_opt(3);
        rho         = param_opt(5);
        kappa_x       = param_opt(6);
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
            WelfCostPrem(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
            WelfCostPrem(k) = NaN;
        end
    end
    
    for k=1:length(rho_grid)
        rho  = rho_grid(k);
        kappa_pie   = param_opt(1);
        kappa_pieW  = param_opt(2);
        kappa_y     = param_opt(3);
        kappa_prem  = param_opt(4);
        kappa_x       = param_opt(6);
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
            WelfCostRho(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
            WelfCostRho(k) = NaN;
        end
    end
    
    
    for k=1:length(kappa_x_grid)
        kappa_x  = kappa_x_grid(k);
        kappa_pie   = param_opt(1);
        kappa_pieW  = param_opt(2);
        kappa_y     = param_opt(3);
        kappa_prem  = param_opt(4);
        rho         = param_opt(5);
        save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
        clear global M_ oo_
        %dynare GK_Nom_CapU noclearall
        try
            GK_Nom_CapU
            W_pos=strmatch('Wf',M_.endo_names,'exact');
            welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
            WelfCostKappx(k) =100*(1-exp((welfare_grid-WelfOutput(1,10))*(1-betta)));
        catch
            WelfCostKappx(k) = NaN;
        end
    end
    
    %     for k=1:length(kappa_x_grid)
    %         kappa_x  = kappa_x_grid(k);
    %         kappa_pie   = 1.5;
    %         kappa_pieW  = 0;
    %         kappa_y     = 0;
    %         kappa_prem  = 0;
    %         rho         = 0;
    %         save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','kappa_x','rho','-append');
    %         clear global M_ oo_
    %         %dynare GK_Nom_CapU noclearall
    %         try
    %             GK_Nom_CapU
    %             W_pos=strmatch('Wf',M_.endo_names,'exact');
    %             welfare_grid=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos));
    %             WelfCostKappx(k) =100*(1-exp((welfare_grid-WelfOutput(1,2))*(1-betta)));
    %         catch
    %             WelfCostKappx(k) = NaN;
    %         end;
    %     end;
    
    oldfolder=cd('Figures\Sticky Wages and Sticky Prices');
    
    %WelfCostY(4) = NaN;
    wfc1=figure('Name','Welfare Cost Coefficients');
    plot(kappa_pie_grid,WelfCostPie,'-','LineWidth',2)
    hold on
    plot(kappa_pieW_grid,WelfCostPieW,'-.','LineWidth',2)
    plot(kappa_pie_grid,WelfCostY,'o-','LineWidth',2,'MarkerSize',8)
    %plot(rho_grid,WelfCostRho,'--','LineWidth',2)
    line([0 5],[0 0]);
    hold off
    ylabel('welfare costs in percent')
    xlabel('coefficient value')
    title('Welfare Costs Deviating from Optimal Values')
    legend('\kappa_\pi','\kappa_{\pi_w}','\kappa_y')
    %legend('WelfCostPie','WelfCostPieW','WelfCostY','WelfCostRho')
    saveas(wfc1,'WelfareCostsCoefficients','fig')
    saveas(wfc1,'WelfareCostsCoefficients','eps')
    
    wfc2=figure('Name','Welfare Cost Indexation');
    plot(rho_grid,WelfCostRho,'LineWidth',2)
    hold on
    line([0 1],[0 0]);
    ylabel('welfare costs in percent')
    xlabel('coefficient value')
    title('Welfare Costs Deviating from Optimal Indexation')
    saveas(wfc2,'WelfareCostsIndex','fig')
    saveas(wfc2,'WelfareCostsIndex','eps')
    
    wfc3=figure('Name','Welfare Cost Weight on Asset Growth');
    plot(kappa_x_grid,WelfCostKappx,'LineWidth',2)
    ylabel('welfare costs in percent')
    xlabel('coefficient value')
    title('Welfare Costs Deviating from Optimal Asset Growth')
    saveas(wfc3,'WelfareCostsX','fig')
    saveas(wfc3,'WelfareCostsX','eps')
    
    
    %     wfc3=figure('Name','Welfare Cost Weight on Asset Growth');
    %     plot(kappa_x_grid,WelfCostKappx,'LineWidth',2); hold on;
    %     line([0 5],[0 0]);
    %     ylabel('welfare costs in percent')
    %     xlabel('coefficient value')
    %     title('Welfare Costs of Asset Growth for \kappa_\pi=1.5')
    %     saveas(wfc3,'WelfareCostsX_2','fig')
    %     saveas(wfc3,'WelfareCostsX_2','eps')
    
    cd(oldfolder);
else
    disp('No robustness checks are printed!')
end

%% RUN_COMPARATIV_STATICS
if comp_static == 1
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Run Comparative Statics!')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%')
    
    
    
    %     dynare GK_Nom_CapU
    %     h=0;
    %     gamma_p=0;
    %     gamma_w=0;
    %     thetaw  = 4;
    %     epsilon = 4;
    %     save('params','h','gamma_p','gamma_w','thetaw','epsilon','-append');
    %     [ys]  = CapU_steadystate;
    %     K_ss        = ys(1);
    %     b           = ys(55);
    %     delta_c     = ys(56);
    %     save('params','K_ss','b','delta_c','irf_plot','-append');
    %     [ys]  = GK_Nom_CapU_steadystate_CS;
    %     Wf_FF       = ys(87);
    %
    %     oldfolder=cd('C:\Users\tobi\ownCloud\PhD\GK_Nom_Wages\Code\GK_NoFinFriction\NK_Cap_pw');
    %     dynare NK_NoFinFric
    %     h=0;
    %     gamma_p=0;
    %     gamma_w=0;
    %     thetaw  = 4;
    %     epsilon = 4;
    %     save('params','h','gamma_p','gamma_w','thetaw','epsilon','-append');
    %     [ys]  = CapU_steadystate;
    %     K_ss        = ys(1);
    %     b           = ys(33);
    %     delta_c     = ys(34);
    %     save('params','K_ss','b','delta_c','irf_plot','-append');
    %     [ys]=NK_NoFinFric_steadystate_CS;
    %     WF_NoFrictions = ys(end);
    %
    %
    %
    %     h=0.8;
    %     gamma_p=0;
    %     gamma_w=0;
    %     thetaw  = 4;
    %     epsilon = 4;
    %     save('params','h','gamma_p','gamma_w','thetaw','epsilon','-append');
    %     [ys]  = CapU_steadystate;
    %     K_ss        = ys(1);
    %     b           = ys(33);
    %     delta_c     = ys(34);
    %     save('params','K_ss','b','delta_c','irf_plot','-append');
    %     [ys]=NK_NoFinFric_steadystate_CS;
    %     WF_Habits = ys(end);
    %
    %     h=0;
    %     gamma_p=0;
    %     gamma_w=0;
    %     thetaw  = 4;
    %     epsilon = 4;
    %     save('params','h','gamma_p','gamma_w','thetaw','epsilon','-append');
    %     [ys]  = CapU_steadystate;
    %     K_ss        = ys(1);
    %     b           = ys(33);
    %     delta_c     = ys(34);
    %     save('params','K_ss','b','delta_c','irf_plot','-append');
    %     [ys]=NK_NoFinFric_steadystate_CS;
    %     WF_MonLabor = ys(end);
    %
    %
    %
    %     h=0;
    %     gamma_p=0;
    %     gamma_w=0;
    %     thetaw  = 4;
    %     epsilon = 4;
    %     save('params','h','gamma_p','gamma_w','thetaw','epsilon','-append');
    %     [ys]  = CapU_steadystate;
    %     K_ss        = ys(1);
    %     b           = ys(33);
    %     delta_c     = ys(34);
    %     save('params','K_ss','b','delta_c','irf_plot','-append');
    %     [ys]=NK_NoFinFric_steadystate_CS;
    %     WF_MonGoods = ys(end);
    %
    %
    %     WF_frictions_vec = [WF_NoFrictions WF_Habits WF_MonLabor WF_MonGoods Wf_FF];
    %
    %     Names = {'No Frictions';'Habits';'Monopolistics Labor';'Monopolisitcs Goods';'Financial Frictions'};
    %
    %
    %     for j=1:size(WF_frictions_vec,2)
    %         WelfCost(j) =100*(1-exp((WF_frictions_vec(1,j)-WF_frictions_vec(1,1))*(1-betta)));
    %     end;
    %
    %     WelfareCosts =WelfCost';
    %     T = table(WF_frictions_vec',WelfareCosts,'RowNames',Names)
    %
    %     cd(oldfolder)
    
    lambda_grid     = 0.1:0.05:0.8;
    Y_vec           = zeros(length(lambda_grid),1);
    MC_vec          = zeros(length(lambda_grid),1);
    K_vec           = zeros(length(lambda_grid),1);
    Premium_vec     = zeros(length(lambda_grid),1);
    MPK_vec         = zeros(length(lambda_grid),1);
    z_vec           = zeros(length(lambda_grid),1);
    Wf_vec          = zeros(length(lambda_grid),1);
    phi_vec         = zeros(length(lambda_grid),1);
    MPL_vec         = zeros(length(lambda_grid),1);
    
    lambda = 0.6;
    dynare GK_Nom_CapU
    for k=1:length(lambda_grid)
        lambda = lambda_grid(k);
        save('params','K_ss','b','delta_c','lambda','irf_plot','-append');
        [ys]  = CapU_steadystate;
        K_ss        = ys(1);
        b           = ys(55);
        delta_c     = ys(56);
        save('params','K_ss','b','delta_c','lambda','irf_plot','-append');
        [ys]  = GK_Nom_CapU_steadystate_CS;
        Y_vec(k)        = ys(1);
        MC_vec(k)       = ys(5);
        K_vec(k)        = ys(22);
        Premium_vec(k)  = ys(35);
        MPK_vec(k)      = ys(56);
        z_vec(k)        = ys(33);
        Wf_vec(k)       = ys(87);
        phi_vec(k)      = ys(32);
        MPL_vec(k)      = ys(55);
    end
    data_vec = [Y_vec  K_vec  Premium_vec  MPK_vec  z_vec  Wf_vec  phi_vec MPL_vec];
    var_names = {'Output','Capital','Premium','Marg. Prod. Capital','Bank Net Worth Growth','Welfare','Leverage','Marg. Prod. Labor'};
    
    spaceH=0.03;
    spaceV=0.05;
    marTop=0.03;
    marBot=0.03;
    marginL=0.1;
    margin=0.1;
    padding=0.0;
    %margin=0.03;
    %marginL=0.075;
    figure('Name','Comparative Statics','NumberTitle','off');
    for j=1:8
        sh=subaxis(2,4,j,'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
        hold on
        plot(lambda_grid,data_vec(:,j),'r-','LineWidth',2)
        plot([0.1 0.8],[data_vec(1,j) data_vec(end,j)],'k-','LineWidth',1)
        axis tight
        hold off
        %         if j >  length(irf.name)-3;
        %             xlabel('quarters','interpreter','latex', 'FontSize', 10)
        %         end
        %         set(gca,'fontsize',6)
        %         if j<10
        %             ylabel('$ \% $dev. from St.St','interpreter','latex', 'FontSize', 10)
        %         else
        %             ylabel('dev. from St.St','interpreter','latex', 'FontSize', 10)
        %         end;
        title(var_names(j),'interpreter','latex', 'FontSize', 12)
        
        if j == 5
            l={'Steady State Value','45 Degree Line'} ;
            
            lh=legend(l,'Orientation','horizontal');
            sp=get(sh,'position');
            set(lh,'position',[sp(1),sp(2)-.12,sp(3),.1]);
        end
        %         if j == length(irf.name)-1
        %             l=legend([h1 h2 h3 h4],'Optimal SR','Benchmark','Pure Inflation Targeting','Wage Targeting and Smoothing','Orientation','horizontal');
        %             set(l,'interpreter','latex','FontSize',14)
        %         end
    end
    
    irf_plot = 1;
    
    Y_std_vec           = zeros(length(lambda_grid),1);
    C_std_vec           = zeros(length(lambda_grid),1);
    L_std_vec           = zeros(length(lambda_grid),1);
    K_std_vec           = zeros(length(lambda_grid),1);
    w_std_vec           = zeros(length(lambda_grid),1);
    phic_std_vec        = zeros(length(lambda_grid),1);
    Wf_std_vec          = zeros(length(lambda_grid),1);
    pie_std_vec         = zeros(length(lambda_grid),1);
    Pie_W_std_vec       = zeros(length(lambda_grid),1);
    Premium_std_vec     = zeros(length(lambda_grid),1);
    
    lambda = 0.381;
    save('params','lambda','irf_plot','-append');
    [ys]  = CapU_steadystate;
    K_ss        = ys(1);
    b           = ys(55);
    delta_c     = ys(56);
    kappa_pie   = 1.5;
    kappa_pieW  = 0;
    kappa_y     = 0.125;
    kappa_prem  = 0;% 1.2
    rho         = 0.8;
    kappa_x       = 0;
    irf_plot = 1;
    save('params','K_ss','b','delta_c','kappa_pie','kappa_x','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','-append')
    save('params','K_ss','b','delta_c','lambda','irf_plot','-append');
    dynare GK_Nom_CapU
    for k=1:length(lambda_grid)
        lambda = lambda_grid(k);
        save('params','K_ss','b','delta_c','lambda','irf_plot','-append');
        [ys]  = CapU_steadystate;
        K_ss        = ys(1);
        b           = ys(55);
        delta_c     = ys(56);
        save('params','K_ss','b','delta_c','lambda','irf_plot','-append');
        GK_Nom_CapU;
        Y_std_vec(k)        = oo_.var(1,1);
        C_std_vec(k)        = oo_.var(2,2);
        L_std_vec(k)        = oo_.var(3,3);
        w_std_vec(k)        = oo_.var(4,4);
        K_std_vec(k)        = oo_.var(5,5);
        phic_std_vec(k)     = oo_.var(6,6);
        pie_std_vec(k)      = oo_.var(7,7);
        Pie_W_std_vec(k)    = oo_.var(8,8);
        Premium_std_vec(k)  = oo_.var(9,9);
        Wf_std_vec(k)       = oo_.var(12,12);
    end
    data_vec = [Y_std_vec  C_std_vec  L_std_vec  w_std_vec  K_std_vec  phic_std_vec  pie_std_vec Pie_W_std_vec Premium_std_vec Wf_std_vec];
    var_names = {'Output','Consumption','Labor','Wages','Capital','Leveragw','Inflation','wage Inflation','Premium','Welfare'};
    
    spaceH=0.03;
    spaceV=0.05;
    marTop=0.03;
    marBot=0.03;
    marginL=0.1;
    margin=0.1;
    padding=0.0;
    %margin=0.03;
    %marginL=0.075;
    figure('Name','Comparative Statics, Standard Deviations','NumberTitle','off');
    for j=1:size(data_vec,2)
        sh=subaxis(2,5,j,'SpacingHoriz', spaceH,'SpacingVert',spaceV, 'PL',padding,'PR',padding,'mt',marTop,'mb',marBot,'ML',marginL,'MR',margin);
        hold on
        plot(lambda_grid,data_vec(:,j),'r-','LineWidth',2)
        %plot([0.1 0.8],[data_vec(j,1) data_vec(j,end)],'k-','LineWidth',1)
        axis tight
        hold off
        %         if j >  length(irf.name)-3;
        %             xlabel('quarters','interpreter','latex', 'FontSize', 10)
        %         end
        %         set(gca,'fontsize',6)
        %         if j<10
        %             ylabel('$ \% $dev. from St.St','interpreter','latex', 'FontSize', 10)
        %         else
        %             ylabel('dev. from St.St','interpreter','latex', 'FontSize', 10)
        %         end;
        title(var_names(j),'interpreter','latex', 'FontSize', 12)
        
        %         if j == 5
        %             l={'Steady State Value','45 Degree Line'} ;
        %
        %             lh=legend(l,'Orientation','horizontal');
        %             sp=get(sh,'position');
        %             set(lh,'position',[sp(1),sp(2)-.12,sp(3),.1]);
        %         end
        %         if j == length(irf.name)-1
        %             l=legend([h1 h2 h3 h4],'Optimal SR','Benchmark','Pure Inflation Targeting','Wage Targeting and Smoothing','Orientation','horizontal');
        %             set(l,'interpreter','latex','FontSize',14)
        %         end
    end
    
else
    disp('No comparative statics are printed')
end


%% RUN Correlations
if run_corr == 1
    %cd(oldfolder);
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Compute Correlation between Variables!')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    lambda_grid     = 0.1:0.05:0.6;
    gamam_p_grid    = 0:0.05:0.85;
    
    corr_Int_Pie       = zeros(length(lambda_grid),length(gamam_p_grid));
    corr_Y_Int       = zeros(length(lambda_grid),length(gamam_p_grid));
    corr_R_Pie       = zeros(length(lambda_grid),length(gamam_p_grid));
    corr_Int_K       = zeros(length(lambda_grid),length(gamam_p_grid));
    corr_R_K       = zeros(length(lambda_grid),length(gamam_p_grid));
    corr_Gap_Inf    = zeros(length(lambda_grid),length(gamam_p_grid));
    corr_Gap_InfW = zeros(length(lambda_grid),length(gamam_p_grid));
    corr_Inf_InfW = zeros(length(lambda_grid),length(gamam_p_grid));
    
    irf_plot = 1;
    kappa_pie   = 0;
    kappa_pieW  = 1.5;
    kappa_y     = 0;
    rho         = 0.75;
    save('params','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','-append');
    
    
    dynare GK_Nom_CapU
    for k = 1:length(lambda_grid)
        lambda = lambda_grid(k);
        for kk = 1:length(gamam_p_grid)
            gamma_p = gamam_p_grid(kk);
            save('params','K_ss','b','delta_c','lambda','gamma_p','irf_plot','-append');
            [ys]  = CapU_steadystate;
            K_ss        = ys(1);
            b           = ys(55);
            delta_c     = ys(56);
            save('params','K_ss','b','delta_c','lambda','gamma_p','irf_plot','-append');
            GK_Nom_CapU
            corr_Int_Pie(k,kk)  = oo_.var(8,9)/(sqrt(oo_.var(8,8))*sqrt(oo_.var(9,9)));
            corr_Y_Int(k,kk)    = oo_.var(8,1)/(sqrt(oo_.var(8,8))*sqrt(oo_.var(1,1)));
            corr_R_Pie(k,kk)    = oo_.var(11,9)/(sqrt(oo_.var(11,11))*sqrt(oo_.var(9,9)));
            corr_Int_K(k,kk)    = oo_.var(8,5)/(sqrt(oo_.var(8,8))*sqrt(oo_.var(5,5)));
            corr_R_K(k,kk)      = oo_.var(5,11)/(sqrt(oo_.var(5,5))*sqrt(oo_.var(11,11)));
            corr_Gap_Inf(k,kk)  = oo_.var(9,16)/(sqrt(oo_.var(9,9))*sqrt(oo_.var(16,16)));
            corr_Gap_InfW(k,kk)  = oo_.var(10,16)/(sqrt(oo_.var(10,10))*sqrt(oo_.var(16,16)));
            corr_Inf_InfW(k,kk)  = oo_.var(10,9)/(sqrt(oo_.var(10,10))*sqrt(oo_.var(9,9)));
        end
    end
    
    c1=figure('Name','Correlation Policy Rate and Inflation','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_Int_Pie');
    rotate3d on
    axis tight
    title('Correlation Policy Rate and Inflation','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    %zlabel('Welfare')
    %az = 0;
    %el = 90;
    %view(az, el);
    c2=figure('Name','Correlation Policy Rate and Output','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_Y_Int');
    rotate3d on
    axis tight
    title('Correlation Policy Rate and Output','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    
    c3=figure('Name','Real Rate and Inflation','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_R_Pie');
    rotate3d on
    axis tight
    title('Correlation Real Rate and Inflation','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    
    c4=figure('Name','Correlation Policy Rate and Capital','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_Int_K');
    rotate3d on
    axis tight
    title('Correlation Policy Rate and Capital','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    
    c5=figure('Name','Correlation Real Rate and Capital','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_R_K');
    rotate3d on
    axis tight
    title('Correlation Real Rate and Capital','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    
    c6=figure('Name','Correlation Output Gap and Inflation','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_Gap_Inf');
    rotate3d on
    axis tight
    title('Correlation Output Gap and Inflation','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    
    
    c7=figure('Name','Correlation Output Gap and Wage Inflation','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_Gap_InfW');
    rotate3d on
    axis tight
    title('Correlation Output Gap and Wage Inflation','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    
    
    c8=figure('Name','Correlation Wage Inflation and Inflation','NumberTitle','off');
    [X,Y] = meshgrid(lambda_grid,gamam_p_grid);
    surf(X,Y,corr_Inf_InfW');
    rotate3d on
    axis tight
    title('Correlation Wage Inflation and Inflation','interpreter','latex', 'FontSize', 12)
    xlabel('\lambda')
    ylabel('\gamma_p')
    
else
    disp('No Correlations are calculated')
end

%% Policy Frontier

if pol_front == 1
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Print Policy Frontiers!')
    disp('%%%%%%%%%%%%%%%%%%%%%%%')
    
    
    irf_plot    = 0;
    osr_ind     = 1;
    kappa_pie   = 1.5;
    kappa_y     = 0.125;
    kappa_pieW  = 0;
    kappa_prem  = 0;
    rho         = 0.8;
    
    
    
    save('params','K_ss','b','delta_c','kappa_pie','kappa_pieW','kappa_y','kappa_prem','rho','irf_plot','osr_ind','-append');
    
    dynare GK_Nom_CapU_OSR_PiPieW noclearall
    
    weight_pie_vec = 0:0.01:1;
    
    frontier_PiPieW     = zeros(length(weight_pie_vec),4);
    frontier_PiY        = zeros(length(weight_pie_vec),4);
    frontier_PiWY       = zeros(length(weight_pie_vec),4);
    
    % Policy frontier for stabilizing price inflation and wage inflation
    for k=1:length(weight_pie_vec)
        weight_pie = weight_pie_vec(k);
        save('params','weight_pie','-append');
        GK_Nom_CapU_OSR_PiPieW
        frontier_PiPieW(k,1)=oo_.var(1,1);% price inflation
        frontier_PiPieW(k,2)=oo_.var(2,2);% outputgap
        frontier_PiPieW(k,3)=oo_.var(3,3);% Wage Inflation
        frontier_PiPieW(k,4)=oo_.var(4,4);% Premium
    end
    
    % Policy frontier for stabilizing price inflation and output gap
    dynare GK_Nom_CapU_OSR_PiY noclearall
    weight_pie_vec = 0:0.01:1;
    for k=1:length(weight_pie_vec)
        weight_pie = weight_pie_vec(k);
        save('params','weight_pie','-append');
        GK_Nom_CapU_OSR_PiY
        frontier_PiY(k,1)=oo_.var(1,1);% price inflation
        frontier_PiY(k,2)=oo_.var(2,2);% outputgap
        frontier_PiY(k,3)=oo_.var(3,3);% Wage Inflation
        frontier_PiY(k,4)=oo_.var(4,4);% Premium
    end
    
    % Policy frontier for stabilizing wage inflation and output gap
    dynare GK_Nom_CapU_OSR_PiWY noclearall
    weight_pie_vec = 0:0.01:1;
    for k=1:length(weight_pie_vec)
        weight_pie = weight_pie_vec(k);
        save('params','weight_pie','-append');
        GK_Nom_CapU_OSR_PiWY
        frontier_PiWY(k,1)=oo_.var(1,1);% price inflation
        frontier_PiWY(k,2)=oo_.var(2,2);% outputgap
        frontier_PiWY(k,3)=oo_.var(3,3);% Wage Inflation
        frontier_PiWY(k,4)=oo_.var(4,4);% Premium
    end
    
    oldfolder=cd('Figures\Sticky Wages and Sticky Prices\Pol_Front_w_kappa_x');
    pf1=figure('Name','Policy Frontier Pie Ygap','NumberTitle','off');
    plot(frontier_PiPieW(1:end,1),frontier_PiPieW(1:end,2),'bo','LineWidth',2,'MarkerSize',5)
    hold on
    plot(frontier_PiY(1:end,1),frontier_PiY(1:end,2),'rx','LineWidth',2,'MarkerSize',5)
    plot(frontier_PiWY(1:end,1),frontier_PiWY(1:end,2),'g+','LineWidth',2,'MarkerSize',5)
    axis tight
    %title('Policy Frontier Benchmark TR','FontSize', 12)
    l=legend('Stab. \Pi \Pi^W','Stab. \Pi Y_{gap}','Stab. \Pi^W Y_{gap}');%,'fontweight','bold','FontSize', 12);
    l.FontWeight = 'bold';
    l.FontSize = 12;
    xlabel('\sigma^2(\Pi)','fontweight','bold','FontSize', 12)
    ylabel('\sigma^2 (Y_{gap})','fontweight','bold','FontSize', 12)
    saveas(pf1,'Front_pie_ygap','fig')
    saveas(pf1,'Front_pie_ygap.eps','epsc')
    
    pf2=figure('Name','Policy Frontier Pie PieW','NumberTitle','off');
    plot(frontier_PiPieW(1:end,1),frontier_PiPieW(1:end,3),'bo','LineWidth',2,'MarkerSize',5)
    hold on
    plot(frontier_PiY(1:end,1),frontier_PiY(1:end,3),'rx','LineWidth',2,'MarkerSize',5)
    plot(frontier_PiWY(1:end,1),frontier_PiWY(1:end,3),'g+','LineWidth',2,'MarkerSize',5)
    axis tight
    %title('Policy Frontier Benchmark TR','FontSize', 12)
l=legend('Stab. \Pi \Pi^W','Stab. \Pi Y_{gap}','Stab. \Pi^W Y_{gap}');%,'fontweight','bold','FontSize', 12);
    l.FontWeight = 'bold';
    l.FontSize = 12;    xlabel('\sigma^2(\Pi)','fontweight','bold','FontSize', 12)
    ylabel('\sigma^2(\Pi^{W})','fontweight','bold','FontSize', 12)
    saveas(pf2,'Front_pie_pieW','fig')
    saveas(pf2,'Front_pie_pieW.eps','epsc')
    
    pf3=figure('Name','Policy Frontier PieW Ygap','NumberTitle','off');
    plot(frontier_PiPieW(1:end,3),frontier_PiPieW(1:end,2),'bo','LineWidth',2,'MarkerSize',5)
    hold on
    plot(frontier_PiY(1:end,3),frontier_PiY(1:end,2),'rx','LineWidth',2,'MarkerSize',5)
    plot(frontier_PiWY(1:end,3),frontier_PiWY(1:end,2),'g+','LineWidth',2,'MarkerSize',5)
    axis tight
    title('Policy Frontier Benchmark TR','FontSize', 12)
l=legend('Stab. \Pi \Pi^W','Stab. \Pi Y_{gap}','Stab. \Pi^W Y_{gap}');%,'fontweight','bold','FontSize', 12);
    l.FontWeight = 'bold';
    l.FontSize = 12;    xlabel('\sigma^2(\Pi^W)','fontweight','bold','FontSize', 12)
    ylabel('\sigma^2(y_{gap})','fontweight','bold','FontSize', 12)
    saveas(pf3,'Front_pieW_ygap','fig')
    saveas(pf3,'Front_pieW_ygap.eps','epsc')
    
    pf4=figure('Name','Policy Frontier PieW Prem','NumberTitle','off');
    plot(frontier_PiPieW(1:end,3),frontier_PiPieW(1:end,4),'bo','LineWidth',2,'MarkerSize',5)
    hold on
    plot(frontier_PiY(1:end,3),frontier_PiY(1:end,4),'rx','LineWidth',2,'MarkerSize',5)
    plot(frontier_PiWY(1:end,3),frontier_PiWY(1:end,4),'g+','LineWidth',2,'MarkerSize',5)
    axis tight
    %title('Policy Frontier Benchmark TR','FontSize', 12)
l=legend('Stab. \Pi \Pi^W','Stab. \Pi Y_{gap}','Stab. \Pi^W Y_{gap}');%,'fontweight','bold','FontSize', 12);
    l.FontWeight = 'bold';
    l.FontSize = 12;    xlabel('\sigma^2(\Pi^W)','fontweight','bold','FontSize', 12)
    ylabel('\sigma^2(Prem)','fontweight','bold','FontSize', 12)
    saveas(pf4,'Front_pieW_prem','fig')
    saveas(pf4,'Front_pieW_prem.eps','epsc')
    
    pf5=figure('Name','Policy Frontier Pie Prem','NumberTitle','off');
    plot(frontier_PiPieW(1:end,1),frontier_PiPieW(1:end,4),'bo','LineWidth',2,'MarkerSize',5)
    hold on
    plot(frontier_PiY(1:end,1),frontier_PiY(1:end,4),'rx','LineWidth',2,'MarkerSize',5)
    plot(frontier_PiWY(1:end,1),frontier_PiWY(1:end,4),'g+','LineWidth',2,'MarkerSize',5)
    axis tight
    %title('Policy Frontier Benchmark TR','FontSize', 12)
l=legend('Stab. \Pi \Pi^W','Stab. \Pi Y_{gap}','Stab. \Pi^W Y_{gap}');%,'fontweight','bold','FontSize', 12);
    l.FontWeight = 'bold';
    l.FontSize = 12;    xlabel('\sigma^2(\Pi)','fontweight','bold','FontSize', 12)
    ylabel('\sigma^2(Prem)','fontweight','bold','FontSize', 12)
    saveas(pf5,'Front_pie_prem','fig')
    saveas(pf5,'Front_pie_prem.eps','epsc')
    
    cd(oldfolder)
else
    disp('No Policy Frontiers are printed')
end
