clc;
clear('all');
rng('default');
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMETERS = {};
PARAMETERS.SYSTEM = 2; % 0-> CSTR; 1-> FITZHUGH-NAGUMO; 2-> ENZYMATIC
PARAMETERS.CREATE_PDF = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (PARAMETERS.SYSTEM == 1) % 1-> FITZHUGH-NAGUMO
    PARAMETERS.SAMPLING_TIME = 1e-2;
    PARAMETERS.TOTAL_TIME = 50;
    PARAMETERS.GLOBAL_GAIN = 1.0;    
    PARAMETERS.ETA = 0.20531;
    PARAMETERS.EPSILON = 0.2966;
    PARAMETERS.BETA = 1.96;
    PARAMETERS.ALFA = 0.50;    
    PARAMETERS.x1 = 1.0656;
    PARAMETERS.x2 = 2.6903;
    PARAMETERS.x3 = fitzhugh_nagumo_model_gV(PARAMETERS.x1, PARAMETERS); 
    PARAMETERS.x1_MAX = 4.0;
    PARAMETERS.x2_MAX = 5.0;
    PARAMETERS.x3_MAX = PARAMETERS.BETA*PARAMETERS.x1_MAX; 
    PARAMETERS.DOT_X1_MAX = PARAMETERS.x1_MAX + (PARAMETERS.x1_MAX^3/3) + PARAMETERS.x2_MAX;
    PARAMETERS.DOT_X2_MAX = PARAMETERS.EPSILON*(PARAMETERS.x3_MAX + PARAMETERS.x2_MAX + PARAMETERS.ETA);
    PARAMETERS.DOT_X3_MAX = PARAMETERS.BETA*PARAMETERS.DOT_X1_MAX; 
    PARAMETERS.DOT_F1_MAX = PARAMETERS.DOT_X2_MAX;
    PARAMETERS.DOT_F2_MAX = PARAMETERS.EPSILON*PARAMETERS.BETA*PARAMETERS.DOT_X3_MAX;     
elseif (PARAMETERS.SYSTEM == 0) % 0-> CSTR
    PARAMETERS.SAMPLING_TIME = 1e-2;
    PARAMETERS.TOTAL_TIME = 100;
    PARAMETERS.GLOBAL_GAIN = 1.0;        
    PARAMETERS.K = 69;
    PARAMETERS.U1 = 1;
    PARAMETERS.U2 = 0;
    PARAMETERS.U3 = -1;
    PARAMETERS.Q = 1;
    PARAMETERS.QC = 0.28;
    PARAMETERS.DELTA1 = 10;
    PARAMETERS.GAMMA = 2.0; 
    PARAMETERS.x1 = 0.12;
    PARAMETERS.x2 = 2.67;
    PARAMETERS.x3 = 0.58;
    PARAMETERS.x1_MAX = 5.0;
    PARAMETERS.x2_MAX = 4.0;
    PARAMETERS.x3_MAX = 4.0; 
    PARAMETERS.DOT_X2_MAX = PARAMETERS.Q*(PARAMETERS.U2+PARAMETERS.x2_MAX) + 0.3*(PARAMETERS.x2_MAX+PARAMETERS.x1_MAX) + 0.576*PARAMETERS.x3_MAX*exp((PARAMETERS.GAMMA*PARAMETERS.x2_MAX)/(PARAMETERS.GAMMA+PARAMETERS.x2_MAX));
    PARAMETERS.DOT_X3_MAX = PARAMETERS.DELTA1*PARAMETERS.QC*(PARAMETERS.U3+PARAMETERS.x1_MAX) + 3*(PARAMETERS.x2_MAX+PARAMETERS.x1_MAX); 
    PARAMETERS.DOT_F1_MAX = 3*PARAMETERS.DOT_X2_MAX;
    PARAMETERS.DOT_F2_MAX = 0.576*(PARAMETERS.DOT_X3_MAX*exp((PARAMETERS.GAMMA*PARAMETERS.x2_MAX)/(PARAMETERS.GAMMA+PARAMETERS.x2_MAX)));
    PARAMETERS.DOT_F3_MAX = 10.0;
 elseif (PARAMETERS.SYSTEM == 2) % 2-> ENZYMATIC
    PARAMETERS.SAMPLING_TIME = 1e-2;
    PARAMETERS.TOTAL_TIME = 50;
    PARAMETERS.GLOBAL_GAIN = 1.0;        
    PARAMETERS.a1 = 360;
    PARAMETERS.k1 = 1;
    PARAMETERS.A1 = 43;
    PARAMETERS.b1 = 1;
    PARAMETERS.alfa1 = 1;
    PARAMETERS.beta1 = 0.6;
    PARAMETERS.gamma1 = 1;
    PARAMETERS.delta1 = 0.8; 
    PARAMETERS.exp1 = 10; 
    PARAMETERS.x1 = 0.0;
    PARAMETERS.x2 = 0.0;
    PARAMETERS.x3 = 0.0;    
    PARAMETERS.x1_MAX = 8.0;
    PARAMETERS.x2_MAX = 8.0;
    PARAMETERS.x3_MAX = 8.0;
    PARAMETERS.DOT_X2_MAX = PARAMETERS.gamma1*PARAMETERS.x3_MAX+PARAMETERS.delta1*PARAMETERS.x2_MAX;
    PARAMETERS.DOT_X3_MAX = PARAMETERS.alfa1*PARAMETERS.x1_MAX+PARAMETERS.beta1*PARAMETERS.x3_MAX;
    PARAMETERS.DOT_F1_MAX = PARAMETERS.a1*((10*PARAMETERS.k1*PARAMETERS.DOT_X2_MAX)/((PARAMETERS.A1)^2));
    PARAMETERS.DOT_F2_MAX = PARAMETERS.gamma1*PARAMETERS.DOT_X3_MAX;  
    PARAMETERS.DOT_F3_MAX = PARAMETERS.DOT_F2_MAX;
    PARAMETERS.alfa_MAX = 5;
end
PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate
simulation_data = run_simulation(PARAMETERS);

% Plot
plot_simulation(simulation_data, PARAMETERS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulation 
function simulation_data = run_simulation(PARAMETERS)
    % Simulation time
    simulation_time = 0:PARAMETERS.SAMPLING_TIME:PARAMETERS.TOTAL_TIME-PARAMETERS.SAMPLING_TIME;
    simulation_steps = size(simulation_time, 2);
    simulation_data = zeros(simulation_steps, 20); % time, x1, x2, u, p, est_x2
        
    % Inital plant states
    x1 = PARAMETERS.x1;
    x2 = PARAMETERS.x2;
    x3 = PARAMETERS.x3; 
    
    % Initial control states
    w_z1 = 0;
    w_z2 = 0.0;
    w_z3 = 0.0;
    diff_est_x1 = x1;
    diff_est_x2 = 0.0;    
    diff_est_x3 = 0.0;
    z_x1 = 0.0;
    z_x2 = 0.0;
    
    % Simulation
    time = 0.0;
    for simulation_step = 1:simulation_steps        
        % Plant
        if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
            % Disturbance value
            if time < 2
                d = 0;
            elseif time < 20
                d = 1;
            elseif time < 30
                d = -(time/10)+3;
            elseif time < 30+8*pi
                d = 2*sin((time-30)/4);
            elseif time < 30+17*pi
                d = 0.7*sin((time-30));
            else
                d = 0;
            end
            
            % Dynamic model
            [x1_new, x2_new, x3_new, f1, g1, f2, g2, f3, g3] = cstr_model(x1, x2, x3, d, PARAMETERS);
            
            % Diferenciation operation (g1->f1)
            [z_x1_new, est_f1, diff_est_x1_new, w_z1_new] = differenciation(x1, g1, diff_est_x1, w_z1, PARAMETERS.SAMPLING_TIME, 1/PARAMETERS.SAMPLING_TIME, PARAMETERS.DOT_F1_MAX);
            est_x1 = diff_est_x1;
            
            % Function inversion (f1->x2->f_2)
            est_x2 = est_f1/3.0;
            if (simulation_step == 1)
                diff_est_x2 = est_x2;
            end
            f2_est = cstr_model_f2(x1, est_x2, PARAMETERS);
            
            % Diferenciation operation (f2->g2)
            [z_x2_new, g2_est, diff_est_x2_new, w_z2_new] = differenciation(est_x2, f2_est, diff_est_x2, w_z2, PARAMETERS.SAMPLING_TIME, 1/PARAMETERS.SAMPLING_TIME, PARAMETERS.DOT_F2_MAX);
            
            % Function inversion (g2->x3)
            est_x3 = g2_est / (0.576*cstr_model_F(est_x2, PARAMETERS));
            est_x3 = min(abs(est_x3), PARAMETERS.x3_MAX)*sign(est_x3);

            % Estimation of external disturbance. Diferenciation operation (g3->f3)            
            [z_x3_new, est_d, diff_est_x3_new, w_z3_new] = differenciation(est_x3, cstr_model_f3(est_x2, est_x3, PARAMETERS), diff_est_x3, w_z3, PARAMETERS.SAMPLING_TIME, 1/PARAMETERS.SAMPLING_TIME, PARAMETERS.DOT_F3_MAX);
            est_d = min(abs(est_d), PARAMETERS.x3_MAX)*sign(est_d);
            
        elseif (PARAMETERS.SYSTEM == 1) % 1-> FITZHUGH-NAGUMO
            
            % Dynamic model
            [x1_new, x2_new, x3_new, f1, g1, f2, g2] = fitzhugh_nagumo_model(x1, x2, x3, PARAMETERS);

            % Diferenciation operation (g1->f1)
            [z_x1_new, est_f1, diff_est_x1_new, w_z1_new] = differenciation(x1, g1, diff_est_x1, w_z1, PARAMETERS.SAMPLING_TIME, 1/PARAMETERS.SAMPLING_TIME, PARAMETERS.DOT_F1_MAX);
            est_x1 = diff_est_x1;
            
            % Function inversion (f1->x2->g_2)
            est_x2 = -est_f1;
            if (simulation_step == 1)
                diff_est_x2 = est_x2;
            end
            g2_est = fitzhugh_nagumo_model_g2(est_x2, PARAMETERS);
            
            % Diferenciation operation (g2->f2)
            [z_x2_new, f2_est, diff_est_x2_new, w_z2_new] = differenciation(est_x2, g2_est, diff_est_x2, w_z2, PARAMETERS.SAMPLING_TIME, 1/PARAMETERS.SAMPLING_TIME, PARAMETERS.DOT_F2_MAX);
            
            % Function inversion (g2->x3)
            est_x3 = f2_est / PARAMETERS.EPSILON;  
            est_x3 = min(abs(est_x3), PARAMETERS.x3_MAX)*sign(est_x3);
            
        elseif (PARAMETERS.SYSTEM == 2) % 2-> ENZYMATIC
            PARAMETERS.alfa1 = 1 + 0.25*exp(-1*cos(3*time) - 1.5*sin(7*time));
            
            % Dynamic model
            [x1_new, x2_new, x3_new, f1, g1, f2, g2] = enzymatic_model(x1, x2, x3, PARAMETERS);
   
            % Diferenciation operation (g1->f1)
            [z_x1_new, est_f1, diff_est_x1_new, w_z1_new] = differenciation(x1, g1, diff_est_x1, w_z1, PARAMETERS.SAMPLING_TIME, 1/PARAMETERS.SAMPLING_TIME, PARAMETERS.DOT_F1_MAX);
            est_x1 = diff_est_x1;
            if sign(est_f1) == 0
                est_f1 = max(abs(est_f1), 1e-8);
            else
                est_f1 = max(abs(est_f1), 1e-8)*sign(est_f1);
            end
            
            % Function inversion (f1->x2->g_2)
            est_x2 = (((PARAMETERS.a1/est_f1)-PARAMETERS.A1)/PARAMETERS.k1)^(1.0/PARAMETERS.exp1);           
            if (simulation_step == 1)
                diff_est_x2 = est_x2;
            end
            
            % Diferenciation operation (g2->f2)
            g2_est = enzymatic_g2(est_x2, PARAMETERS);
            [z_x2_new, f2_est, diff_est_x2_new, w_z2_new] = differenciation(est_x2, g2_est, diff_est_x2, w_z2, PARAMETERS.SAMPLING_TIME, 1 /PARAMETERS.SAMPLING_TIME,  PARAMETERS.DOT_F2_MAX);
                        
            % Function inversion (g2->x3)
            est_x3 = f2_est / PARAMETERS.gamma1;  
            est_x3 = min(abs(est_x3), PARAMETERS.x3_MAX)*sign(est_x3);

            % Estimation of inknown parameter alfa 
            g3_est = enzymatic_g3(est_x3, PARAMETERS);
            [z_x3_new, est_f3, diff_est_x3_new, w_z3_new] = differenciation(est_x3, g3_est, diff_est_x3, w_z3, PARAMETERS.SAMPLING_TIME, 1/PARAMETERS.SAMPLING_TIME, PARAMETERS.DOT_F3_MAX);
            est_alfa = est_f3/x1; 
            est_alfa = min(abs(est_alfa), PARAMETERS.alfa_MAX)*sign(est_alfa);            
            est_x2_condition = PARAMETERS.A1+(PARAMETERS.k1*x2^PARAMETERS.exp1);
     
        end

        % Save data
        simulation_data(simulation_step, 1) = time;
        simulation_data(simulation_step, 2) = x1;
        simulation_data(simulation_step, 3) = x2;
        simulation_data(simulation_step, 4) = x3;  
        simulation_data(simulation_step, 5) = est_x1;
        simulation_data(simulation_step, 6) = est_x2;
        simulation_data(simulation_step, 7) = est_x3;
        if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
            simulation_data(simulation_step, 8) = d;
            simulation_data(simulation_step, 9) = est_d;
        elseif (PARAMETERS.SYSTEM == 2) % 0-> ENZYMATIC
            simulation_data(simulation_step, 8) = est_alfa;
            simulation_data(simulation_step, 9) = PARAMETERS.alfa1;
            simulation_data(simulation_step, 10) = est_x2_condition;
        end             
        simulation_data(simulation_step, 15) = z_x1;
        simulation_data(simulation_step, 16) = z_x2;        
        
        % Update data
        x1 = x1_new;
        x2 = x2_new;
        x3 = x3_new;
        z_x1 = z_x1_new;
        z_x2 = z_x2_new;
        diff_est_x1 = diff_est_x1_new;
        diff_est_x2 = diff_est_x2_new;
        if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
            diff_est_x3 = diff_est_x3_new;
        elseif (PARAMETERS.SYSTEM == 2) % 0-> ENZYMATIC
            diff_est_x3 = diff_est_x3_new;
        end
        w_z1 = w_z1_new;        
        w_z2 = w_z2_new;
        if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
            w_z3 = w_z3_new;
        elseif (PARAMETERS.SYSTEM == 2) % 0-> ENZYMATIC
            w_z3 = w_z3_new;
        end
        time = time + PARAMETERS.SAMPLING_TIME;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIFFERENCTIATOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differenciation function dx = f + d -> estimation of d
function [z, est_d, est_x_new, w_z_new] = differenciation(x, f, est_x, w_z, tau, lambda, beta)
    nu = (((tau/pi)^2)*(beta/2))^(1/3);
    gamma = beta / nu;
    z = x - est_x;    
    est_d = lambda*z + w_z;
    dot_est_x = f + est_d;
    dot_w_z = (((lambda^2))/4.0)*z + saturation(gamma*z, beta);
    est_x_new = est_x + dot_est_x*tau;
    w_z_new = w_z + dot_w_z*tau;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F] = cstr_model_F(x2, PARAMETERS)
    F = exp((PARAMETERS.GAMMA*x2)/(PARAMETERS.GAMMA+x2));
end
function [f2] = cstr_model_f2(x1, x2, PARAMETERS)
    f2 = PARAMETERS.Q*(PARAMETERS.U2-x2) - 0.3*(x2-x1);
end
function [f3] = cstr_model_f3(x2, x3, PARAMETERS)
    f3 = PARAMETERS.Q*(PARAMETERS.U1-x3) - 0.072*cstr_model_F(x2, PARAMETERS)*x3;
end

function [x1_new, x2_new, x3_new, f1, g1, f2, g2, f3, g3] = cstr_model(x1, x2, x3, d, PARAMETERS)
    f1 = 3*x2;
    g1 = PARAMETERS.DELTA1*PARAMETERS.QC*(PARAMETERS.U3-x1) - 3*x1;
    f2 = cstr_model_f2(x1, x2, PARAMETERS);
    g2 = 0.576*cstr_model_F(x2, PARAMETERS)*x3;   
    f3 = cstr_model_f3(x2, x3, PARAMETERS);
    g3 = d;    
    dot_x1 = f1 + g1;    
    dot_x2 = f2 + g2;
    dot_x3 = f3 + g3;
    x1_new = x1 + dot_x1*PARAMETERS.SAMPLING_TIME;
    x2_new = x2 + dot_x2*PARAMETERS.SAMPLING_TIME;
    x3_new = x3 + dot_x3*PARAMETERS.SAMPLING_TIME;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENZYMATIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g3] = enzymatic_g3(x3, PARAMETERS)
    g3 = -PARAMETERS.beta1*x3;
end

function [g2] = enzymatic_g2(x2, PARAMETERS)
    g2 = -PARAMETERS.delta1*x2;
end

function [x1_new, x2_new, x3_new, f1, g1, f2, g2] = enzymatic_model(x1, x2, x3, PARAMETERS)
    f1 = PARAMETERS.a1/(PARAMETERS.A1+(PARAMETERS.k1*x2^PARAMETERS.exp1));
    g1 = -PARAMETERS.b1*x1;  
    g2 = enzymatic_g2(x2, PARAMETERS); 
    f2 = PARAMETERS.gamma1*x3; 
    g3 = enzymatic_g3(x3, PARAMETERS);
    f3 = PARAMETERS.alfa1*x1;    
    dot_x1 = f1 + g1;    
    dot_x2 = f2 + g2;
    dot_x3 = f3 + g3;   
    x1_new = x1 + dot_x1*PARAMETERS.SAMPLING_TIME;
    x2_new = x2 + dot_x2*PARAMETERS.SAMPLING_TIME;
    x3_new = x3 + dot_x3*PARAMETERS.SAMPLING_TIME;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FITZHUGH-NAGUMO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gV] = fitzhugh_nagumo_model_gV(x1, PARAMETERS)
    if (x1 > 0)
        gV = PARAMETERS.BETA*x1;
    else
        gV = PARAMETERS.ALFA*x1;
    end
end

function [g2] = fitzhugh_nagumo_model_g2(x2, PARAMETERS)
    g2 = -PARAMETERS.EPSILON*(x2 + PARAMETERS.ETA);
end

function [x1_new, x2_new, x3_new, f1, g1, f2, g2] = fitzhugh_nagumo_model(x1, x2, x3, PARAMETERS)
    f1 = -x2;
    g1 = x1 - ((x1^3)/3);  
    g2 = fitzhugh_nagumo_model_g2(x2, PARAMETERS); 
    f2 = PARAMETERS.EPSILON*x3;    
    dot_x1 = f1 + g1;    
    dot_x2 = f2 + g2;
    x1_new = x1 + dot_x1*PARAMETERS.SAMPLING_TIME;
    x2_new = x2 + dot_x2*PARAMETERS.SAMPLING_TIME;
    x3_new = fitzhugh_nagumo_model_gV(x1, PARAMETERS);
end

% Plot simulation 
function  plot_simulation(simulation_data, PARAMETERS)
    plot_font_size = 12;
    if (PARAMETERS.SYSTEM == 0) % 1-> CSTR
        figure(1);
        clf(1);
        subplot(3,1,1);
        plot(simulation_data(:,1), simulation_data(:,3) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,6) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$T_r(t)$ [K]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $T_r(t)$','FontSize', plot_font_size,'Interpreter','latex');
        legend('$T_r(t)$','$\hat{T}_r(t)$', 'Interpreter','latex','FontSize', plot_font_size);        
        subplot(3,1,2);
        plot(simulation_data(:,1), simulation_data(:,4) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,7) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$C_r(t)$ [mol/l]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $C_r(t)$','FontSize', plot_font_size,'Interpreter','latex');
        legend('$C_r(t)$','$\hat{C}_r(t)$', 'Interpreter','latex','FontSize', plot_font_size);        
        subplot(3,1,3);
        plot(simulation_data(:,1), simulation_data(:,8) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,9) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$d(t)$ [mol/l]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $d(t)$','FontSize', plot_font_size,'Interpreter','latex');
        legend('$d(t)$','$\hat{d}(t)$', 'Interpreter','latex','FontSize', plot_font_size);
        if PARAMETERS.CREATE_PDF
            export_fig('../MANUSCRIPT/GRAPHICS/cstr_states.pdf', '-transparent');
        end
        figure(2);
        clf(2);
        T_END = 3.0;
        subplot(3,1,1);
        plot(simulation_data(:,1), simulation_data(:,3) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,6) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$T_r(t)$ [K]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Detail for $t \in[0, 3.0]$ of $\hat{T}_r(t)$','FontSize', plot_font_size,'Interpreter','latex');
        xlim([0.0, T_END]); 
        legend('$T_r(t)$','$\hat{T}_r(t)$', 'Interpreter','latex','FontSize', plot_font_size);
        subplot(3,1,2);
        plot(simulation_data(:,1), simulation_data(:,4) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,7) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$C_r(t)$ [mol/l]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Detail for $t \in[0, 3.0]$ of $\hat{C}_r(t)$','FontSize', plot_font_size,'Interpreter','latex');
        xlim([0.0, T_END]); 
        legend('$C_r(t)$','$\hat{C}_r(t)$', 'Interpreter','latex','FontSize', plot_font_size);
        subplot(3,1,3);
        plot(simulation_data(:,1), simulation_data(:,8) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,9) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$d(t)$ [mol/l]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Detail for $t \in[0, 3.0]$ of $\hat{d}(t)$','FontSize', plot_font_size,'Interpreter','latex');
        xlim([0.0, T_END]);
        legend('$d(t)$','$\hat{d}(t)$', 'Interpreter','latex','FontSize', plot_font_size);
        if PARAMETERS.CREATE_PDF
            export_fig('../MANUSCRIPT/GRAPHICS/cstr_states_detail.pdf', '-transparent');
        end
    elseif (PARAMETERS.SYSTEM == 1) % 1-> FITZHUGH-NAGUMO  
        figure(1);
        clf(1);
        subplot(2,1,1);
        plot(simulation_data(:,1), simulation_data(:,3) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,6) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$W(t)$ vs $\hat{W}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $W(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
        legend('W(t)','$\hat{W}(t)$', 'Interpreter','latex','FontSize', plot_font_size);        
        subplot(2,1,2);
        plot(simulation_data(:,1), simulation_data(:,4) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,7) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('g(V) vs $\hat{g}(V)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of g(V)','FontSize', plot_font_size,'Interpreter','latex');
        legend('g(V)','$\hat{g}(V)$', 'Interpreter','latex','FontSize', plot_font_size);
        if PARAMETERS.CREATE_PDF
            export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_states.pdf', '-transparent');
        end 
        figure(2);
        clf(2);
        subplot(2,1,1);
        plot(simulation_data(:,1), simulation_data(:,3) ,'-b', 'LineWidth',1.5);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,6) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$W(t)$ vs $\hat{W}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Detail for $t \in[0, 2.0]$ of $\hat{W}(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
        legend('W(t)','$\hat{W}(t)$', 'Interpreter','latex','FontSize', plot_font_size);
        xlim([0.0, 2.0]);
        subplot(2,1,2);
        plot(simulation_data(:,1), simulation_data(:,4) ,'-b', 'LineWidth',1.5);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,7) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('g(V) vs $\hat{g}(V)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Detail for $t \in[0, 2.0]$ of $\hat{g}(V)$','FontSize', plot_font_size,'Interpreter','latex');   
        xlim([0.0, 2.0]);
        legend('g(V)','$\hat{g}(V)$', 'Interpreter','latex','FontSize', plot_font_size);
        if PARAMETERS.CREATE_PDF
            export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_states_detail.pdf', '-transparent');
        end 
        figure(3);
        clf(3);
        plot(simulation_data(:,2), simulation_data(:,4) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,2), simulation_data(:,7) ,'--r', 'LineWidth',1.0);
        xlabel('V(t)', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('g(V) vs $\hat{g}(V)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Trajectory in plane (V,g(V))','FontSize', plot_font_size,'Interpreter','latex'); 
        legend('(V,g(V))','(V,$\hat{g}(V))$', 'Interpreter','latex', 'Location','northwest','FontSize', plot_font_size);
        if PARAMETERS.CREATE_PDF
            export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_phase_plane.pdf', '-transparent');
        end 
     elseif (PARAMETERS.SYSTEM == 2) % 2-> ENZYMATIC 
        figure(1);
        clf(1);
        subplot(3,1,1);
        plot(simulation_data(:,1), simulation_data(:,3) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,6) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$Z(t)$ vs $\hat{Z}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $Z(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
        legend('Z(t)','$\hat{Z}(t)$', 'Interpreter','latex','FontSize', plot_font_size);        
        subplot(3,1,2);
        plot(simulation_data(:,1), simulation_data(:,4) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,7) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$Y(t)$ vs $\hat{Y}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $Y(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
        legend('Y(t)','$\hat{Y}(t)$', 'Interpreter','latex','FontSize', plot_font_size', 'Location','southeast'); 
        subplot(3,1,3);
        plot(simulation_data(:,1), simulation_data(:,9) ,'-b', 'LineWidth',1.0);
        grid on;
        hold on;
        plot(simulation_data(:,1), simulation_data(:,8) ,'--r', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$\alpha(t)$ vs $\hat{\alpha}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of parameter $\alpha(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
        legend('$\alpha(t)$','$\hat{\alpha}(t)$', 'Interpreter','latex','FontSize', plot_font_size, 'Location','southeast'); 
        if PARAMETERS.CREATE_PDF
            export_fig('../MANUSCRIPT/GRAPHICS/enzymatic.pdf', '-transparent');
        end
    end
end

% Saturation function
function y = saturation(x, x_max)
    if abs(x)> x_max
        y=x_max*sign(x);
    else 
        y=x;
    end
end

