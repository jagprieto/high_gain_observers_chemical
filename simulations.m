clc;
clear('all');
rng('default');
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMETERS = {};
PARAMETERS.CREATE_PDF = 0; 
PARAMETERS.SYSTEM = 0; % 0-> CSTR; 1-> FITZHUGH-NAGUMO
PARAMETERS.NOISE_ACTIVATION = 1;
PARAMETERS.SAMPLING_TIME = 1e-2;
if PARAMETERS.NOISE_ACTIVATION  == 0
    PARAMETERS.NOISE_MODULE_DB = 0;
elseif PARAMETERS.NOISE_ACTIVATION  == 1
    PARAMETERS.NOISE_MODULE_DB = 30;
elseif PARAMETERS.NOISE_ACTIVATION  == 2
    PARAMETERS.NOISE_AMPLITUDE = 0.25;
    PARAMETERS.NOISE_FREQUENCY = 80;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (PARAMETERS.SYSTEM == 1) % 1-> FITZHUGH-NAGUMO
    PARAMETERS.TOTAL_TIME = 50;
    PARAMETERS.GLOBAL_GAIN = 1.0;    
    PARAMETERS.ETA = 0.20531;
    PARAMETERS.EPSILON = 0.2966;
    PARAMETERS.BETA = 1.96;
    PARAMETERS.ALFA = 0.50;    
    PARAMETERS.x1 = 1.0656;
    PARAMETERS.x2 = 2.6903;
    PARAMETERS.x3 = fitzhugh_nagumo_gV(PARAMETERS.x1, PARAMETERS); 
    PARAMETERS.x1_MAX = 4.0;
    PARAMETERS.x2_MAX = 5.0;
    PARAMETERS.x3_MAX = PARAMETERS.BETA*PARAMETERS.x1_MAX; 
elseif (PARAMETERS.SYSTEM == 0) % 0-> CSTR
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
    PARAMETERS.x2_MAX = 5.0;
    PARAMETERS.x3_MAX = 4.0; 
    PARAMETERS.x4_MAX = 5.0; 
end

PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate
simulation_data = run_simulation(PARAMETERS);

% Plot
plot_simulation(simulation_data, PARAMETERS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIFFERENCTIATOR FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differenciation function dx = f + d
function [e, est_d, est_x_new, w_e_new] = differenciation_snsta(x, f, est_x, w_e, omega_c, PARAMETERS)
     e = x - est_x;
     lambda_gain = 2*omega_c;
     gamma_gain = 2/PARAMETERS.SAMPLING_TIME;
     beta_gain = 1*lambda_gain/(gamma_gain);
     dot_w_e = ((lambda_gain^2/(4)))*e + (lambda_gain/2)*beta_gain*gamma_gain*((sech(gamma_gain*e))^2)*e;
     %alfa_gain = lambda_gain + beta_gain*gamma_gain*((sech(gamma_gain*e))^2);
     %dot_w_e = ((alfa_gain^2/(4)))*e + (lambda_gain/4)*beta_gain*tanh(gamma_gain*e);
     w_e_new = w_e + dot_w_e*PARAMETERS.SAMPLING_TIME;          
     dot_est_x = f + lambda_gain*e + beta_gain*tanh(gamma_gain*e) + w_e_new;
     est_x_new = est_x + dot_est_x*PARAMETERS.SAMPLING_TIME;     
     est_d = w_e_new;
end

function [e, est_d, est_x_new, w_e_new, dot_w_e] = differenciation_st(x, f, est_x, w_e, omega_c, PARAMETERS)
     e = x - est_x;
     lambda_1 = 0.5*omega_c;     
     lambda_2 = 2*lambda_1;
     dot_w_e = lambda_2*sign(e);
     w_e_new = w_e + dot_w_e*PARAMETERS.SAMPLING_TIME;
     dot_est_x = f + lambda_1*((abs(e))^0.5)*sign(e) + w_e;     
     est_x_new = est_x + dot_est_x*PARAMETERS.SAMPLING_TIME;     
     est_d = dot_est_x - f; %w_e_new;%
end

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
    
    % Initial SNTA control states
    w_z1 = 0;
    w_z2 = 0.0;
    w_z3 = 0.0;
    diff_est_x1 = 0.0;
    diff_est_x2 = 0.0;    
    diff_est_x3 = 0.0;
    z_x1 = 0.0;
    z_x2 = 0.0;
    
    % Initial ST control states
    st_w_z1 = 0;
    st_w_z2 = 0.0;
    st_w_z3 = 0.0;
    st_diff_est_x1 = 0.0;
    st_diff_est_x2 = 0.0;    
    st_diff_est_x3 = 0.0;
    st_z_x1 = 0.0;
    st_z_x2 = 0.0;
    
    % Cascade observer data
    cascade_observer_state = zeros(4,1);
    cascade_estimated_state = zeros(4,1);
    
    
    % Simulation
    time = 0.0;
    for simulation_step = 1:simulation_steps 
        
        % Read output
        if PARAMETERS.NOISE_ACTIVATION == 1
            y = awgn(x1 , PARAMETERS.NOISE_MODULE_DB, "measured");   
            noise = y - x1;
        elseif PARAMETERS.NOISE_ACTIVATION == 2
            noise = PARAMETERS.NOISE_AMPLITUDE*cos(PARAMETERS.NOISE_FREQUENCY*time);
            noise = noise + 0.83*PARAMETERS.NOISE_AMPLITUDE*sin(1.29*PARAMETERS.NOISE_FREQUENCY*time - 0.14);
            noise = noise + 0.23*PARAMETERS.NOISE_AMPLITUDE*cos(5.12*PARAMETERS.NOISE_FREQUENCY*time + 0.26);          
            noise = noise + 0.65*PARAMETERS.NOISE_AMPLITUDE*sin(3.37*PARAMETERS.NOISE_FREQUENCY*time + 0.36)*exp(cos(1.21*PARAMETERS.NOISE_FREQUENCY*time + 0.13)); 
            y = x1 + noise;
        else
            y = x1;
            noise = 0;
        end
        
        if simulation_step < 2
            diff_est_x1 = y;
            st_diff_est_x1 = y;
            cascade_observer_state(1) = y;
        end
        
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
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNSTA DIFFERENCIATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Crossover frequencies
            convex_gain = 0.80;            
            if PARAMETERS.NOISE_ACTIVATION == 0
                omega_c1 = 0.25/PARAMETERS.SAMPLING_TIME;
                st_gains = 2*[0.25, 0.05, 0.025];
                exponent_gains = [.25, .1];
            elseif PARAMETERS.NOISE_ACTIVATION == 1
                omega_c1 = min(0.25/PARAMETERS.SAMPLING_TIME, 10);
                st_gains = [2*0.5, 0.2, 0.2];
                exponent_gains = [.25, .1];
            elseif PARAMETERS.NOISE_ACTIVATION == 2
                omega_c1 = min(0.25/PARAMETERS.SAMPLING_TIME, 8);
                st_gains = 2*[0.5, 0.25, 0.1];
                exponent_gains = [.25, .1];
            end            
            omega_c2 = omega_c1*((1-convex_gain) + convex_gain*(1.0-exp(-exponent_gains(1)*time)));
            omega_c3 = omega_c1*((1-convex_gain) + convex_gain*(1.0-exp(-exponent_gains(2)*time)));
            
            % Diferenciation operation (g1->f1)
            [z_x1_new, est_f1, diff_est_x1_new, w_z1_new] = differenciation_snsta(y, g1, diff_est_x1, w_z1, omega_c1, PARAMETERS);
            est_x1 = diff_est_x1;
            est_x1 = min(abs(est_x1), PARAMETERS.x2_MAX)*sign(est_x1);

            % Function inversion (f1->x2->f_2)
            est_x2 = cstr_inverse_f1(est_f1,PARAMETERS);%est_f1/3.0;
            est_x2 = min(abs(est_x2), PARAMETERS.x2_MAX)*sign(est_x2);
            if (simulation_step == 1)
                diff_est_x2 = est_x2;
            end
            g2_est = cstr_g2(est_x1, est_x2, PARAMETERS);
            
            % Diferenciation operation (g2->f2)
            [z_x2_new, est_f2, diff_est_x2_new, w_z2_new] = differenciation_snsta(est_x2, g2_est, diff_est_x2, w_z2, omega_c2, PARAMETERS);
            
            % Function inversion (f2->x3)
            est_x3 = cstr_inverse_f2(est_f2, est_x2, PARAMETERS);%est_f2 / (0.576*cstr_model_F(est_x2, PARAMETERS));
            est_x3 = min(abs(est_x3), PARAMETERS.x3_MAX)*sign(est_x3);

            % Estimation of external disturbance. Diferenciation operation (g3->f3)  
            [z_x3_new, est_d, diff_est_x3_new, w_z3_new] = differenciation_snsta(est_x3, cstr_f3(est_x2, est_x3, PARAMETERS), diff_est_x3, w_z3, omega_c3, PARAMETERS);
            est_d = min(abs(est_d), PARAMETERS.x3_MAX)*sign(est_d);
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ST DIFFERENCIATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Diferenciation operation (g1->f1)
            [st_z_x1_new, st_est_f1, st_diff_est_x1_new, st_w_z1_new] = differenciation_st(y, g1, st_diff_est_x1, st_w_z1, st_gains(1)*omega_c1, PARAMETERS);
            st_est_x1 = st_diff_est_x1;
            
            % Function inversion (f1->x2->f_2)
            st_est_x2 = cstr_inverse_f1(st_est_f1,PARAMETERS);
            if (simulation_step == 1)
                st_diff_est_x2 = st_est_x2;
            end
            st_g2_est = cstr_g2(y, st_est_x2, PARAMETERS);
            
            % Diferenciation operation (g2->f2)
            [st_z_x2_new, st_est_f2, st_diff_est_x2_new, st_w_z2_new] = differenciation_st(st_est_x2, st_g2_est, st_diff_est_x2, st_w_z2, st_gains(2)*omega_c2, PARAMETERS);
            
            % Function inversion (f2->x3)
            st_est_x3 = cstr_inverse_f2(st_est_f2, st_est_x2, PARAMETERS);
            st_est_x3 = min(abs(st_est_x3), PARAMETERS.x3_MAX)*sign(st_est_x3);

            % Estimation of external disturbance. Diferenciation operation (g3->f3)  
            [st_z_x3_new, st_est_d, st_diff_est_x3_new, st_w_z3_new] = differenciation_st(st_est_x3, cstr_f3(st_est_x2, st_est_x3, PARAMETERS), st_diff_est_x3, st_w_z3, st_gains(3)*omega_c3, PARAMETERS);
            st_est_d = min(abs(st_est_d), PARAMETERS.x3_MAX)*sign(st_est_d);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH GAIN CASCADE OBSERVER (KHALIL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [cascade_estimated_state, cascade_observer_state] = khalil_cascade_observer(y, cascade_observer_state, PARAMETERS);
            
        elseif (PARAMETERS.SYSTEM == 1) % 1-> FITZHUGH-NAGUMO
            % Crossover frequencies
            convex_gain = 0.80;
            if PARAMETERS.NOISE_ACTIVATION == 0
                exponent_gain = 0.25;
                omega_c1 = 0.25/PARAMETERS.SAMPLING_TIME;
                st_gains = [0.25, 0.1];
                omega_c2 = 1*omega_c1*((1-convex_gain) + convex_gain*(1.0-exp(-exponent_gain*time)));
            elseif PARAMETERS.NOISE_ACTIVATION == 1
                exponent_gain = 0.25;
                omega_c1 = min(0.25/PARAMETERS.SAMPLING_TIME, 8);
                st_gains = [0.25, 0.1];
                omega_c2 = omega_c1*((1-convex_gain) + convex_gain*(1.0-exp(-exponent_gain*time)));
            elseif PARAMETERS.NOISE_ACTIVATION == 2
                exponent_gain = 0.25;
                omega_c1 = min(0.25/PARAMETERS.SAMPLING_TIME, 7);
                st_gains = [0.25, 0.1];
                omega_c2 = omega_c1*((1-convex_gain) + convex_gain*(1.0-exp(-exponent_gain*time)));
            end              
            
  
            % Dynamic model
            [x1_new, x2_new, x3_new, f1, g1, f2, g2] = fitzhugh_nagumo_model(x1, x2, x3, PARAMETERS);
            
            % Diferenciation operation (g1->f1)
            [z_x1_new, est_f1, diff_est_x1_new, w_z1_new] = differenciation_snsta(y, g1, diff_est_x1, w_z1, omega_c1, PARAMETERS);
            est_x1 = diff_est_x1;   
 
            
            % Function inversion (f1->x2->g_2)
            est_x2 = -est_f1;
            if (simulation_step == 1)
                diff_est_x2 = est_x2;
            end
            g2_est = fitzhugh_nagumo_g2(est_x2, PARAMETERS);
            
            % Diferenciation operation (g2->f2)
            [z_x2_new, est_f2, diff_est_x2_new, w_z2_new] = differenciation_snsta(est_x2, g2_est, diff_est_x2, w_z2, omega_c2, PARAMETERS);
            
            % Function inversion (g2->x3)
            est_x3 = est_f2 / PARAMETERS.EPSILON;  
            est_x3 = min(abs(est_x3), PARAMETERS.x3_MAX)*sign(est_x3);
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ST DIFFERENCIATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Diferenciation operation (g1->f1)
            [st_z_x1_new, st_est_f1, st_diff_est_x1_new, st_w_z1_new] = differenciation_st(y, g1, st_diff_est_x1, st_w_z1, st_gains(1)*omega_c1, PARAMETERS);
            st_est_x1 = st_diff_est_x1;
                        
            % Function inversion (f1->x2->g_2)
            st_est_x2 = -st_est_f1;
            if (simulation_step == 1)
                st_diff_est_x2 = st_est_x2;
            end
            st_g2_est = fitzhugh_nagumo_g2(st_est_x2, PARAMETERS);
            
            % Diferenciation operation (g2->f2)
            [st_z_x2_new, st_est_f2, st_diff_est_x2_new, st_w_z2_new] = differenciation_st(st_est_x2, st_g2_est, st_diff_est_x2, st_w_z2, st_gains(2)*omega_c2, PARAMETERS);
         
            % Function inversion (g2->x3)
            st_est_x3 = st_est_f2 / PARAMETERS.EPSILON;  
            st_est_x3 = min(abs(st_est_x3), PARAMETERS.x3_MAX)*sign(st_est_x3);        
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH GAIN CASCADE OBSERVER (KHALIL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [cascade_estimated_state, cascade_observer_state] = khalil_cascade_observer(y, cascade_observer_state, PARAMETERS);
                     
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
        
        simulation_data(simulation_step, 17) = cascade_estimated_state(1);        
        simulation_data(simulation_step, 18) = cascade_estimated_state(2);        
        simulation_data(simulation_step, 19) = cascade_estimated_state(3);        
        simulation_data(simulation_step, 20) = cascade_estimated_state(4);        
        
        simulation_data(simulation_step, 21) = st_est_x1;
        simulation_data(simulation_step, 22) = st_est_x2;
        simulation_data(simulation_step, 23) = st_est_x3;
        if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
            simulation_data(simulation_step, 24) = st_est_d;
        elseif (PARAMETERS.SYSTEM == 2) % 0-> ENZYMATIC
            simulation_data(simulation_step, 24) = st_est_alfa;
            simulation_data(simulation_step, 25) = PARAMETERS.alfa1;
            simulation_data(simulation_step, 26) = st_est_x2_condition;
        end        
        simulation_data(simulation_step, 27) = st_z_x1;
        simulation_data(simulation_step, 28) = st_z_x2;  
        simulation_data(simulation_step, 29) = y;  
        
        % Update state
        x1 = x1_new;
        x2 = x2_new;
        x3 = x3_new;
        
        % Update SNSTA state
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
        
        % Update ST state
        st_z_x1 = st_z_x1_new;
        st_z_x2 = st_z_x2_new;
        st_diff_est_x1 = st_diff_est_x1_new;
        st_diff_est_x2 = st_diff_est_x2_new;
        if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
            st_diff_est_x3 = st_diff_est_x3_new;
        elseif (PARAMETERS.SYSTEM == 2) % 0-> ENZYMATIC
            st_diff_est_x3 = st_diff_est_x3_new;
        end
        st_w_z1 = st_w_z1_new;        
        st_w_z2 = st_w_z2_new;
        if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
            st_w_z3 = st_w_z3_new;
        elseif (PARAMETERS.SYSTEM == 2) % 0-> ENZYMATIC
            st_w_z3 = st_w_z3_new;
        end
        
        % Update time
        time = time + PARAMETERS.SAMPLING_TIME;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH GAIN CASCADE OBSERVER (KHALIL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estimated_system_state, observer_state_new] = khalil_cascade_observer(y, observer_state, PARAMETERS)
 
 scale_factor = 0.015;
 if PARAMETERS.NOISE_ACTIVATION == 1
    epsilon = 10*scale_factor;
 elseif PARAMETERS.NOISE_ACTIVATION == 2
    epsilon = 20*scale_factor;
 else
    epsilon = 1*scale_factor;
 end
  z1 = observer_state(1);
  z2 = observer_state(2);
  z3 = observer_state(3);
  gain = 1; 
  dot_z1 = (1/epsilon)*(z2+2*(y-z1));
  dot_z2 = (1/epsilon)*(y-z1);      
  M2 = gain*PARAMETERS.x2_MAX;
  z2_est = M2*saturation(z2/(M2*epsilon), 1);              
  dot_z3 = (-1/epsilon)*(z3+z2_est);
  M3 = gain*PARAMETERS.x3_MAX;
  z3_est = M3*saturation((z3+z2_est)/(M3*epsilon), 1); 
      
  if (PARAMETERS.SYSTEM == 0) % 0-> CSTR
      
    x1_est = z1;
    g1 = cstr_g1(y, PARAMETERS);
    est_f1 = z2_est - g1;
    x2_est = cstr_inverse_f1(est_f1, PARAMETERS);
    x2_est = saturation(x2_est, PARAMETERS.x2_MAX);      
    g2 = cstr_g2(x1_est, x2_est, PARAMETERS);
    est_f2 = (1/3)*(z3_est - 3*g2 + (PARAMETERS.DELTA1*PARAMETERS.QC+3)*(est_f1+g1));
    x3_est = cstr_inverse_f2(est_f2, x2_est, PARAMETERS);      
    x3_est = saturation(x3_est, PARAMETERS.x3_MAX);
      
      
  elseif (PARAMETERS.SYSTEM == 1) % 1-> FITZHUGH-NAGUMO

    x1_est = z1;
    g1 = fitzhugh_nagumo_g1(y, PARAMETERS);
    est_f1 = z2_est - g1;
    x2_est = fitzhugh_nagumo_inverse_f1(est_f1, PARAMETERS);
    x2_est = saturation(x2_est, PARAMETERS.x2_MAX);      
    g2 = fitzhugh_nagumo_g2(x2_est, PARAMETERS);
    est_f2 =  - z3_est - g2 + z2_est*(1-y^2);
    x3_est = fitzhugh_nagumo_inverse_f2(est_f2, PARAMETERS);      
    x3_est = saturation(x3_est, PARAMETERS.x3_MAX);

  end  
  
  z1 = z1 + dot_z1*PARAMETERS.SAMPLING_TIME;
  z2 = z2 + dot_z2*PARAMETERS.SAMPLING_TIME;
  z3 = z3 + dot_z3*PARAMETERS.SAMPLING_TIME;
  observer_state_new = observer_state;
  observer_state_new(1) = z1;
  observer_state_new(2) = z2;
  observer_state_new(3) = z3;
  estimated_system_state = zeros(4,1);
  estimated_system_state(1) = x1_est;
  estimated_system_state(2) = x2_est; 
  estimated_system_state(3) = x3_est; 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F] = cstr_model_F(x2, PARAMETERS)
    F = exp((PARAMETERS.GAMMA*x2)/(PARAMETERS.GAMMA+x2));
end
function [g2] = cstr_g2(x1, x2, PARAMETERS)
    g2 = PARAMETERS.Q*(PARAMETERS.U2-x2) - 0.3*(x2-x1);
end
function [f2] = cstr_f2(x2, x3, PARAMETERS)
    f2 =  0.576*cstr_model_F(x2, PARAMETERS)*x3;
end

function [f3] = cstr_f3(x2, x3, PARAMETERS)
    f3 = PARAMETERS.Q*(PARAMETERS.U1-x3) - 0.072*cstr_model_F(x2, PARAMETERS)*x3;
end
function [g1] = cstr_g1(x1, PARAMETERS)
    g1 = PARAMETERS.DELTA1*PARAMETERS.QC*(PARAMETERS.U3-x1) - 3*x1;
end
function [f1] = cstr_f1(x2,PARAMETERS)
    f1 = 3*x2;
end

function [x2] = cstr_inverse_f1(f1,PARAMETERS)
    x2 = f1/3;
end
function [x3] = cstr_inverse_f2(f2, x2, PARAMETERS)
    x3 = f2/(0.576*cstr_model_F(x2, PARAMETERS));
end

function [x1_new, x2_new, x3_new, f1, g1, f2, g2, f3, g3] = cstr_model(x1, x2, x3, d, PARAMETERS)
    f1 = cstr_f1(x2,PARAMETERS);
    g1 = cstr_g1(x1, PARAMETERS);
    g2 = cstr_g2(x1, x2, PARAMETERS);
    f2 = cstr_f2(x2, x3, PARAMETERS);  
    f3 = cstr_f3(x2, x3, PARAMETERS);
    g3 = d;    
    dot_x1 = f1 + g1;    
    dot_x2 = f2 + g2;
    dot_x3 = f3 + g3;
    x1_new = x1 + dot_x1*PARAMETERS.SAMPLING_TIME;
    x2_new = x2 + dot_x2*PARAMETERS.SAMPLING_TIME;
    x3_new = x3 + dot_x3*PARAMETERS.SAMPLING_TIME;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FITZHUGH-NAGUMO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gV] = fitzhugh_nagumo_gV(x1, PARAMETERS)
    if (x1 > 0)
        gV = PARAMETERS.BETA*x1;
    else
        gV = PARAMETERS.ALFA*x1;
    end
end

function [g2] = fitzhugh_nagumo_g2(x2, PARAMETERS)
    g2 = -PARAMETERS.EPSILON*(x2 + PARAMETERS.ETA);
end

function [f2] = fitzhugh_nagumo_f2(x3, PARAMETERS)
    f2 = PARAMETERS.EPSILON*x3;
end

function [x3] = fitzhugh_nagumo_inverse_f2(f2, PARAMETERS)
    x3 = f2/PARAMETERS.EPSILON;
end

function [f1] = fitzhugh_nagumo_f1(x2, PARAMETERS)
    f1 = -x2;
end

function [x2] = fitzhugh_nagumo_inverse_f1(f1, PARAMETERS)
    x2 = -f1;
end

function [g1] = fitzhugh_nagumo_g1(x1, PARAMETERS)
    g1 = x1 - ((x1^3)/3);
end

function [x1_new, x2_new, x3_new, f1, g1, f2, g2] = fitzhugh_nagumo_model(x1, x2, x3, PARAMETERS)
    f1 = fitzhugh_nagumo_f1(x2, PARAMETERS);
    g1 = fitzhugh_nagumo_g1(x1, PARAMETERS);  
    g2 = fitzhugh_nagumo_g2(x2, PARAMETERS); 
    f2 = fitzhugh_nagumo_f2(x3, PARAMETERS);    
    dot_x1 = f1 + g1;    
    dot_x2 = f2 + g2;
    x1_new = x1 + dot_x1*PARAMETERS.SAMPLING_TIME;
    x2_new = x2 + dot_x2*PARAMETERS.SAMPLING_TIME;
    x3_new = fitzhugh_nagumo_gV(x1, PARAMETERS);
end



% Saturation function
function y = saturation(x, x_max)
    if abs(x)> x_max
        y=x_max*sign(x);
    else 
        y=x;
    end
end