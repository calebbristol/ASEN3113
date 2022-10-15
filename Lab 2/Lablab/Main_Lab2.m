%%%%%%%%%%%%%%%%%%%%%%%
%  ASEN 3113: Lab 2
%  Group Members:
%  -Caleb Bristol 
%  -Jared Seefried 
%  -Henri Wessels
%  -Justin Pedersen
%  -Sophia Trissel
%
%  Date: 1 April, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%


%% Clear Workspace
clc
clear
close all;


%% Read in Data
Al_a = load('Exp. Lab 2 - Heat Conduction/Aluminum_25V_240mA');
Al_b = load('Exp. Lab 2 - Heat Conduction/Aluminum_28V_269mA');
Br_c = load('Exp. Lab 2 - Heat Conduction/Brass_26V_245mA');
Br_d = load('Exp. Lab 2 - Heat Conduction/Brass_29V_273mA');
St_e = load('Exp. Lab 2 - Heat Conduction/Steel_21V_192mA');


%% Question 1
%
%

    %% Variables
    dist = 1.375:0.5:4.875;

    %% Calculate Slope
    a = polyfit(dist,Al_a(end,3:10),1);
    b = polyfit(dist,Al_b(end,3:10),1);
    c = polyfit(dist,Br_c(end,3:10),1);
    d = polyfit(dist,Br_d(end,3:10),1);
    e = polyfit(dist,St_e(end,3:10),1);

    T_0_a = polyval(a,0);
    T_0_b = polyval(b,0);
    T_0_c = polyval(c,0);
    T_0_d = polyval(d,0);
    T_0_e = polyval(e,0);

    H_a = a(1);
    H_b = b(1);
    H_c = c(1);
    H_d = d(1);
    H_e = e(1);

%% Question 2
%
%


%% Question 3
%
%


%% Question 4
%
%

    %% Establish problem constraints
    alpha_Al = 4.836e-05 * 1550; %[in^2/s]
    alpha_Br = 3.565e-05 * 1550; %[in^2/s]
    alpha_St = 4.05e-06 * 1550; %[in^2/s]

    L = 5; %[in]
    x_ = 4.875; %[in]
    alpha_Al_ = alpha_Al/100:alpha_Al/100:1.1*alpha_Al; %[in^2/s]
    alpha_Br_ = alpha_Br/100:alpha_Br/100:1.1*alpha_Br; %[in^2/s]
    alpha_St_ = alpha_St/100:alpha_St/100:1.1*alpha_St; %[in^2/s]
    t_a = 1:10:10*length(Al_a);
    t_b = 1:10:10*length(Al_b);
    t_c = 1:10:10*length(Br_c);
    t_d = 1:10:10*length(Br_d);
    t_e = 1:10:10*length(St_e);
    
    %% Run Simulation
    %
    % 
    
        %% Aluminum
        temp_an_a = zeros(length(t_a),length(alpha_Al_));
        temp_an_b = zeros(length(t_b),length(alpha_Al_));

        for i = 1:length(t_a)
            for j = 1:length(alpha_Al_)
                temp_an_a(i,j) = u(x_,t_a(i),T_0_a,H_a,alpha_Al_(j),L,1);
            end
        end
        
        for i = 1:length(t_b)
            for j = 1:length(alpha_Al_)
                temp_an_b(i,j) = u(x_,t_b(i),T_0_b,H_b,alpha_Al_(j),L,1);
            end
        end
    
        %% Brass
        temp_an_c = zeros(length(t_c),length(alpha_Br_));
        temp_an_d = zeros(length(t_d),length(alpha_Br_));
        
        for i = 1:length(t_c)
            for j = 1:length(alpha_Br_)
                temp_an_c(i,j) = u(x_,t_c(i),T_0_c,H_c,alpha_Br_(j),L,1);
            end
        end
        
        for i = 1:length(t_d)
            for j = 1:length(alpha_Br_)
                temp_an_d(i,j) = u(x_,t_d(i),T_0_d,H_d,alpha_Br_(j),L,1);
            end
        end
        
        %% Steel
        temp_an_e = zeros(length(t_e),length(alpha_St_));
        
        for i = 1:length(t_e)
            for j = 1:length(alpha_St_)
                temp_an_e(i,j) = u(x_,t_e(i),T_0_e,H_e,alpha_St_(j),L,1);
            end
        end
    
    %% Compare to Experimental Data
    %
    %
    
        %% Aluminum
        temp_dif_a = temp_an_a - repmat(Al_a(:,10),1,length(temp_an_a(1,:)));
        temp_dif_b = temp_an_b - repmat(Al_b(:,10),1,length(temp_an_b(1,:)));
        
        [~,I_a] = min(var(temp_dif_a));
        [~,I_b] = min(var(temp_dif_b));
        
        alpha_a = alpha_Al_(I_a);
        alpha_b = alpha_Al_(I_b);
        
        %% Brass
        temp_dif_c = temp_an_c - repmat(Br_c(:,10),1,length(temp_an_c(1,:)));
        temp_dif_d = temp_an_d - repmat(Br_d(:,10),1,length(temp_an_d(1,:)));
        
        [~,I_c] = min(var(temp_dif_c));
        [~,I_d] = min(var(temp_dif_d));
        
        alpha_c = alpha_Br_(I_c);
        alpha_d = alpha_Br_(I_d);
        
        %% Steel
        temp_dif_e = temp_an_e - repmat(St_e(:,10),1,length(temp_an_e(1,:)));
        
        [~,I_e] = min(var(temp_dif_e));
        
        alpha_e = alpha_St_(I_e);
    
    
    %% Plotting
    figure()
    plot(Al_a(:,1),Al_a(:,10),'r','LineWidth',2); hold on
    for i = 1:length(alpha_Al_)
        plot(t_a,temp_an_a(:,i),'k')
    end
    plot(Al_a(:,1),Al_a(:,10),'r','LineWidth',2);
    title('Temperature vs. Time for Various \alpha Values')
    xlabel('Time [s]')
    ylabel('Predicted Temperature [C]')  
    legend('Experimental','Simulated')
    grid on; grid minor;
    hold off
    
    figure()
    plot(Al_b(:,1),Al_b(:,10),'r','LineWidth',2); hold on
    for i = 1:length(alpha_Al_)
        plot(t_b,temp_an_b(:,i),'k')
    end   
    plot(Al_b(:,1),Al_b(:,10),'r','LineWidth',2);
    title('Temperature vs. Time for Various \alpha Values')
    xlabel('Time [s]')
    ylabel('Predicted Temperature [C]') 
    legend('Experimental','Simulated')
    grid on; grid minor;
    hold off
    
    figure()
    plot(Br_c(:,1),Br_c(:,10),'r','LineWidth',2); hold on
    for i = 1:length(alpha_Al_)
        plot(t_c,temp_an_c(:,i),'k')
    end
    plot(Br_c(:,1),Br_c(:,10),'r','LineWidth',2);
    title('Temperature vs. Time for Various \alpha Values')
    xlabel('Time [s]')
    ylabel('Predicted Temperature [C]')  
    legend('Experimental','Simulated')
    grid on; grid minor;
    hold off
    
    figure()
    plot(Br_d(:,1),Br_d(:,10),'r','LineWidth',2); hold on
    for i = 1:length(alpha_Al_)
        plot(t_d,temp_an_d(:,i),'k')
    end
    plot(Br_d(:,1),Br_d(:,10),'r','LineWidth',2);
    title('Temperature vs. Time for Various \alpha Values')
    xlabel('Time [s]')
    ylabel('Predicted Temperature [C]')  
    legend('Experimental','Simulated')
    grid on; grid minor;
    hold off
    
    figure()
    plot(St_e(:,1),St_e(:,10),'r','LineWidth',2); hold on
    for i = 1:length(alpha_Al_)
        plot(t_e,temp_an_e(:,i),'k')
    end
    plot(St_e(:,1),St_e(:,10),'r','LineWidth',2);
    title('Temperature vs. Time for Various \alpha Values')
    xlabel('Time [s]')
    ylabel('Predicted Temperature [C]')  
    legend('Experimental','Simulated')
    grid on; grid minor;
    hold off


%% Question 5
%
%

    %% Solve For Steady State
    temp_f_a = Al_a(end,10);
    temp_f_b = Al_b(end,10);
    temp_f_c = Br_c(end,10);
    temp_f_d = Br_d(end,10);
    
    t_ss_a = 0;
    t_ss_b = 0;
    t_ss_c = 0;
    t_ss_d = 0;
    
    i = 1;
    while t_ss_a == 0
        if abs((Al_a(i,10)-temp_f_a)/temp_f_a) <= 0.05
            t_ss_a = Al_a(i,1);
        end
        i = i + 1;
    end
    
    i = 1;
    while t_ss_b == 0
        if abs((Al_b(i,10)-temp_f_b)/temp_f_b) <= 0.05
            t_ss_b = Al_b(i,1);
        end
        i = i + 1;
    end
    
    i = 1;
    while t_ss_c == 0
        if abs((Br_c(i,10)-temp_f_c)/temp_f_c) <= 0.05
            t_ss_c = Br_c(i,1);
        end
        i = i + 1;
    end
    
    i = 1;
    while t_ss_d == 0
        if abs((Br_d(i,10)-temp_f_d)/temp_f_d) <= 0.05
            t_ss_d = Br_d(i,1);
        end
        i = i + 1;
    end



%% Functions

function uxt = u(x,t,T_0,H,alpha,L,N)
%% u(x,t)
% Function to solve u(x,t) for a given N
% Used to solve problem 5

    %% Create n dependent variable functions
        b_n = @(n) (8*H*L)/(pi^2*(2*n-1)^2)*(-1)^(n);
        lambda_n = @(n) (2*n-1)*pi/(2*L);

    %% Run Summation
    sum = 0;
    for i = 1:N
        sum = sum + b_n(i)*sin(lambda_n(i)*x)*exp(-(lambda_n(i)^2*alpha*t));
    end
    
    %% Calculate u
    uxt = T_0 + H*x + sum;
end
