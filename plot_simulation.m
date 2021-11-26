
% Plot simulation 
function  plot_simulation(simulation_data, PARAMETERS)
    plot_font_size = 11;
    legend_font_size = 9;
        
    if (PARAMETERS.SYSTEM == 0) % 1-> CSTR  
        figure(1);
        clf(1);
        subplot(4,1,1);
        if PARAMETERS.NOISE_ACTIVATION > 0
            plot(simulation_data(:,1), simulation_data(:,29) ,'--y', 'LineWidth',0.5);
            grid on;
            hold on;
            plot(simulation_data(:,1), simulation_data(:,17) ,'-k', 'LineWidth',1.0);
        else
            plot(simulation_data(:,1), simulation_data(:,17) ,'-k', 'LineWidth',1.0);
            grid on;
            hold on;
        end     
        if PARAMETERS.NOISE_ACTIVATION == 2
            legend_texts = {'$T_c(t)$ (noisy)','KC','SNSTA','$T_c(t)$'};
        elseif PARAMETERS.NOISE_ACTIVATION == 1
            plot(simulation_data(:,1), simulation_data(:,21) ,'-b', 'LineWidth',0.5);
            legend_texts = {'$T_c(t)$ (noisy)','KC','ST','SNSTA','$T_c(t)$'};
        else
            plot(simulation_data(:,1), simulation_data(:,21) ,'-b', 'LineWidth',0.5);
            legend_texts = {'KC','ST','SNSTA','$T_c(t)$'};
        end 
        plot(simulation_data(:,1), simulation_data(:,5) ,'-r', 'LineWidth',1.0);
        plot(simulation_data(:,1), simulation_data(:,2) ,'--c', 'LineWidth',1.0);              
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$T_C(t)$ [K]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $T_C(t)$','FontSize', plot_font_size,'Interpreter','latex');
%         legend(legend_texts, 'Interpreter','latex','FontSize', legend_font_size, 'Location','NorthEast'); 
        
        subplot(4,1,2);
        plot(simulation_data(:,1), simulation_data(:,18) ,'-k', 'LineWidth',1.0);
        grid on;
        hold on;
        if PARAMETERS.NOISE_ACTIVATION < 2
            plot(simulation_data(:,1), simulation_data(:,22) ,'-b', 'LineWidth',0.5);
        end
        plot(simulation_data(:,1), simulation_data(:,6) ,'-r', 'LineWidth',1.0); 
        plot(simulation_data(:,1), simulation_data(:,3) ,'--c', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$T_r(t)$ [K]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $T_r(t)$','FontSize', plot_font_size,'Interpreter','latex');
              
        subplot(4,1,3);
        plot(simulation_data(:,1), simulation_data(:,19) ,'-k', 'LineWidth',1.0);
        grid on;
        hold on;
        if PARAMETERS.NOISE_ACTIVATION < 2
            plot(simulation_data(:,1), simulation_data(:,23) ,'-b', 'LineWidth',0.5);
        end
        plot(simulation_data(:,1), simulation_data(:,7) ,'-r', 'LineWidth',1.0);
        plot(simulation_data(:,1), simulation_data(:,4) ,'--c', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$C_r(t)$ [mol/l]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $C_r(t)$','FontSize', plot_font_size,'Interpreter','latex');    
        
        subplot(4,1,4);
        if PARAMETERS.NOISE_ACTIVATION < 2
            plot(simulation_data(:,1), simulation_data(:,24) ,'-b', 'LineWidth',0.5); 
            grid on;
            hold on;
            plot(simulation_data(:,1), simulation_data(:,9) ,'-r', 'LineWidth',1.0); 
        else
            plot(simulation_data(:,1), simulation_data(:,9) ,'-r', 'LineWidth',1.0); 
            grid on;
            hold on;
        end        
        plot(simulation_data(:,1), simulation_data(:,8) ,'--c', 'LineWidth',1.05);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$d(t)$ [mol/l]', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $d(t)$','FontSize', plot_font_size,'Interpreter','latex');
        
        if PARAMETERS.CREATE_PDF
            if PARAMETERS.NOISE_ACTIVATION == 0
                export_fig('../MANUSCRIPT/GRAPHICS/cstr_states_no_noise.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 1
                export_fig('../MANUSCRIPT/GRAPHICS/cstr_states_gaussian_noise.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 2
                export_fig('../MANUSCRIPT/GRAPHICS/cstr_states_harmonic_noise.pdf', '-transparent', '-nocrop');
            end
        end
        
%         
%         T_END = 2.0;
        if PARAMETERS.CREATE_PDF
            if PARAMETERS.NOISE_ACTIVATION == 0
                export_fig('../MANUSCRIPT/GRAPHICS/cstr_states_no_noise_detail.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 1
                export_fig('../MANUSCRIPT/GRAPHICS/cstr_states_gaussian_noise_detail.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 2
                export_fig('../MANUSCRIPT/GRAPHICS/cstr_states_harmonic_noise_detail.pdf', '-transparent', '-nocrop');
            end
        end
        
    elseif (PARAMETERS.SYSTEM == 1) % 1-> FITZHUGH-NAGUMO  
        figure(1);
        clf(1);
        subplot(3,1,1);
        if PARAMETERS.NOISE_ACTIVATION > 0
            plot(simulation_data(:,1), simulation_data(:,29) ,'--y', 'LineWidth',0.5);
            grid on;
            hold on;
            plot(simulation_data(:,1), simulation_data(:,17) ,'-k', 'LineWidth',1.0);
        else
            plot(simulation_data(:,1), simulation_data(:,17) ,'-k', 'LineWidth',1.0);
            grid on;
            hold on;
        end     
        if PARAMETERS.NOISE_ACTIVATION == 2
            legend_texts = {'$V(t)$ (noisy)','KC','SNSTA','$V(t)$'};
        elseif PARAMETERS.NOISE_ACTIVATION == 1
            plot(simulation_data(:,1), simulation_data(:,21) ,'-b', 'LineWidth',0.5);
            legend_texts = {'$V(t)$ (noisy)','KC','ST','SNSTA','$V(t)$'};
        else
            
            plot(simulation_data(:,1), simulation_data(:,21) ,'-b', 'LineWidth',0.5);
            legend_texts = {'KC','ST','SNSTA','$V(t)$'};
        end 
        plot(simulation_data(:,1), simulation_data(:,5) ,'-r', 'LineWidth',1.0);
        plot(simulation_data(:,1), simulation_data(:,2) ,'--c', 'LineWidth',1.0); 
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$V(t)$ vs $\hat{V}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $V(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
%         legend(legend_texts, 'Interpreter','latex','FontSize', legend_font_size, 'Location','NorthEast'); 
        
        subplot(3,1,2);
        plot(simulation_data(:,1), simulation_data(:,18) ,'-k', 'LineWidth',1.0);
        grid on;
        hold on;
        if PARAMETERS.NOISE_ACTIVATION < 2
            plot(simulation_data(:,1), simulation_data(:,22) ,'-b', 'LineWidth',0.5);
        end
        plot(simulation_data(:,1), simulation_data(:,6) ,'-r', 'LineWidth',1.0); 
        plot(simulation_data(:,1), simulation_data(:,3) ,'--c', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$W(t)$ vs $\hat{W}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $W(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
        
        subplot(3,1,3);
        plot(simulation_data(:,1), simulation_data(:,19) ,'-k', 'LineWidth',1.0);
        grid on;
        hold on;
        if PARAMETERS.NOISE_ACTIVATION < 2
            plot(simulation_data(:,1), simulation_data(:,23) ,'-b', 'LineWidth',0.5);
        end
        plot(simulation_data(:,1), simulation_data(:,7) ,'-r', 'LineWidth',1.0);
        plot(simulation_data(:,1), simulation_data(:,4) ,'--c', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('g(V) vs $\hat{g}(V)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of g(V)','FontSize', plot_font_size,'Interpreter','latex');
        
        if PARAMETERS.CREATE_PDF
            if PARAMETERS.NOISE_ACTIVATION == 0
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_no_noise.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 1
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_gaussian_noise.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 2
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_harmonic_noise.pdf', '-transparent', '-nocrop');
            end
        end
        
        T_END = 2.0;
        figure(2);
        clf(2);
        subplot(3,1,1);
        if PARAMETERS.NOISE_ACTIVATION > 0
            plot(simulation_data(:,1), simulation_data(:,29) ,'--y', 'LineWidth',0.5);
            grid on;
            hold on;
            plot(simulation_data(:,1), simulation_data(:,17) ,'-k', 'LineWidth',1.0);
        else
            plot(simulation_data(:,1), simulation_data(:,17) ,'-k', 'LineWidth',1.0);
            grid on;
            hold on;
        end     
        if PARAMETERS.NOISE_ACTIVATION == 2
            legend_texts = {'$V(t)$ (noisy)','KC','SNSTA','$V(t)$'};
        elseif PARAMETERS.NOISE_ACTIVATION == 1
            plot(simulation_data(:,1), simulation_data(:,21) ,'-b', 'LineWidth',0.5);
            legend_texts = {'$V(t)$ (noisy)','KC','ST','SNSTA','$V(t)$'};
        else
            
            plot(simulation_data(:,1), simulation_data(:,21) ,'-b', 'LineWidth',0.5);
            legend_texts = {'KC','ST','SNSTA','$V(t)$'};
        end 
        plot(simulation_data(:,1), simulation_data(:,5) ,'-r', 'LineWidth',1.0);
        plot(simulation_data(:,1), simulation_data(:,2) ,'--c', 'LineWidth',1.0); 
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$V(t)$ vs $\hat{V}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $V(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
%         legend(legend_texts, 'Interpreter','latex','FontSize', legend_font_size, 'Location','NorthEast'); 
        xlim([0.0, T_END]);
        subplot(3,1,2);
        plot(simulation_data(:,1), simulation_data(:,18) ,'-k', 'LineWidth',1.0);
        grid on;
        hold on;
        if PARAMETERS.NOISE_ACTIVATION < 2
            plot(simulation_data(:,1), simulation_data(:,22) ,'-b', 'LineWidth',0.5);
        end
        plot(simulation_data(:,1), simulation_data(:,6) ,'-r', 'LineWidth',1.0); 
        plot(simulation_data(:,1), simulation_data(:,3) ,'--c', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('$W(t)$ vs $\hat{W}(t)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of $W(t)$','FontSize', plot_font_size,'Interpreter','latex'); 
        xlim([0.0, T_END]);
        subplot(3,1,3);
        plot(simulation_data(:,1), simulation_data(:,19) ,'-k', 'LineWidth',1.0);
        grid on;
        hold on;
        if PARAMETERS.NOISE_ACTIVATION < 2
            plot(simulation_data(:,1), simulation_data(:,23) ,'-b', 'LineWidth',0.5);
        end
        plot(simulation_data(:,1), simulation_data(:,7) ,'-r', 'LineWidth',1.0);
        plot(simulation_data(:,1), simulation_data(:,4) ,'--c', 'LineWidth',1.0);
        xlabel('Time [s]', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('g(V) vs $\hat{g}(V)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Estimation of g(V)','FontSize', plot_font_size,'Interpreter','latex');
        xlim([0.0, T_END]);

        if PARAMETERS.CREATE_PDF
            if PARAMETERS.NOISE_ACTIVATION == 0
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_no_noise_detail.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 1
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_gaussian_noise_detail.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 2
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_harmonic_noise_detail.pdf', '-transparent', '-nocrop');
            end
        end

        figure(3);
        clf(3);
        plot(simulation_data(:,2), simulation_data(:,4) ,'-c', 'LineWidth',1.0);
        grid on;
        hold on;
        if PARAMETERS.NOISE_ACTIVATION  < 2
            plot(simulation_data(:,2), simulation_data(:,23) ,'-b', 'LineWidth',0.5);
        end
        plot(simulation_data(:,2), simulation_data(:,19) ,'-k', 'LineWidth',1.0);
        plot(simulation_data(:,2), simulation_data(:,7) ,'-r', 'LineWidth',1.0);
        xlabel('V(t)', 'FontSize', plot_font_size,'Interpreter','latex');
        ylabel('g(V) vs $\hat{g}(V)$', 'FontSize', plot_font_size,'Interpreter','latex');
        title('Trajectory in plane (V,g(V))','FontSize', plot_font_size,'Interpreter','latex'); 
        if PARAMETERS.CREATE_PDF
            if PARAMETERS.NOISE_ACTIVATION == 0
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_no_noise_phase_plane.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 1
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_gaussian_noise_phase_plane.pdf', '-transparent', '-nocrop');
            elseif PARAMETERS.NOISE_ACTIVATION == 2
                export_fig('../MANUSCRIPT/GRAPHICS/fitzhugh_nagumo_harmonic_noise_phase_plane.pdf', '-transparent', '-nocrop');
            end
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

