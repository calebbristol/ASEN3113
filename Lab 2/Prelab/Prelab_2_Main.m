%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ASEN 3113 Prelab 2
%  Author: Caleb Bristol
%  Date: 11 March, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clean Workspace
clear
clc
close all;


%% Problem 3
  

    %% Establish Measurements
    T = [18.53 22.47 26.87 30.05 35.87 38.56 41.50 46.26];
    dist = 1.375:0.5:4.875;

    %% Best Fit
    [a] = polyfit(dist,T,1);

    T_0 = polyval(a,0);
    T_pred = polyval(a,dist);

    H = a(2);
    
    %% Compare Results
    figure()
    scatter(dist,T); hold on
    scatter(dist,T_pred)
    title('Predicted vs. Measured Temperature Values')
    xlabel('Distance from Last Thermocouple [in]')
    ylabel('Temperature [C]')
    legend('Measured','Predicted')
    grid on; grid minor;
    hold off


%% Problem 5

    
    %% Define Problem Constraints
    % H = problem 3 [K/in]
    alpha = 4.836e-05; %[m^2/s]
    alpha = alpha * 1550; %[in^2/s]
    L = 5; %[in]
    x_ = 4.875; %[in]
    N = 0:10;
    
    %% Run Simulation
    for i = 1:length(N)
        u_t1(i) = u(x_,1,T_0,H,alpha,L,N(i));
        u_t1000(i) = u(x_,1000,T_0,H,alpha,L,N(i));
    end
    
    %% Plot Results
    figure()
    plot(N,u_t1); hold on
    plot(N,u_t1000)
    title('Predicted Temperature Based on Number of Summations')
    xlabel('Number of Summations [N]')
    ylabel('Predicted Temperature [C]')
    legend('t = 1 [s]','t = 1000 [s]')    
    grid on; grid minor;
    hold off
    
    %% Determine Fourier Number
    F_0_1 = alpha*1/L^2;
    F_0_1000 = alpha*1000/L^2;
    
    
%% Problem 6


    %% Establish problem constraints
    alpha_ = [alpha/2 alpha 2*alpha]; %[in^2/s]
    t = 1:1000;
    
    %% Run Simulation
    %
    % Yeah, this is gonna be a big one
    temp = zeros(length(t),length(alpha_));
    
    for i = 1:length(t)
        for j = 1:length(alpha_)
            temp(i,j) = u(x_,t(i),T_0,H,alpha_(j),L,1);
        end
    end
    
    %% Plotting
    figure()
    for i = 1:length(alpha_)
        plot(t,temp(:,i)); hold on
    end
    title('Temperature vs. Time for Various \alpha Values')
    xlabel('Time [s]')
    ylabel('Predicted Temperature [C]')  
    legend('\alpha = \alpha/2','\alpha = \alpha','\alpha = 2\alpha')
    grid on; grid minor;
    hold off
    
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
    