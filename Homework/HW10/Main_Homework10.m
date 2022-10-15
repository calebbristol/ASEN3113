%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ASEN 3113 Homework 10
%  Author: Caleb Bristol
%  Date: 27 April, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all;
clc


%% Helpful Functions / Constants
%
%

    %% Constants
    sigma = 5.67e-08;
    g = 9.81; %[m/s^2]

    %% Grashof Number
    Gr = @(g,beta,T_s,T_inf,L_c,v) (g*beta*(T_s-T_inf)*L_c^3)/v^2;
    
    %% Nusselt Vertical Plate
    Nu_VP = @(Ra,Pr) (0.825 + (0.387*Ra^(1/6)) / (1 + (0.492/Pr)^(9/16))^(8/27))^2;


%% Problem 20-16
%
%

    %% Given
    L = 0.2; %[m]
    t = 25e-03; %[m]
    k_p = 15; %[W/mK]
    T_s_h = 100 + 273; %[K]
    T_c = 7 + 273; %[K]
    T_guess = 12.5 + 273; %[K]
    
    % Sourced from internet because the textbook wasn't helpful
    v = 1e-03^2; %[m^2/s]
    beta = 210e-06; %[1/K]
    alpha = 1.4558e-07; %[m^2/s]
    
    
    %% Solve for Grashof Number
    Gr_1 = Gr(g,beta,T_guess,T_c,L,v);
    
    %% Solve for Prandtl Number
    Pr_1 = v/alpha;
    
    %% Solve for Rayleigh Number
    Ra_1 = Gr_1 * Pr_1;
    
    %% Solve for Nusselt Number
    Nu_1 = Nu_VP(Ra_1,Pr_1);
    
    %% Compare to Condunction Heat Transfer
    h = Nu_1 * k_p / L;
    
    q_conv = h*(T_guess - T_c);
    
    q_cond = k_p/t*(T_s_h - T_guess);
    
    % After iterating between a couple different values, the final
    % temperature settled on was as follows:
    fprintf('Problem 20-16: \n \n')
    fprintf('T_s [deg C]: \n')
    disp(T_guess-273)


%% Problem 20-30
%
%

    %% Given
    A_pcb = 0.15 * 0.20; %[m^2]
    T_inf = 20; %[deg C]
    Q_dis = 8; %[W]
    epsilon = 0.8;
    
    %% Properties of Surrounding Air (Table A-22)
    %
    % Here we are going to take an initial guess that  T_s = 45 [deg C]
    % T_av = 32.5 [deg C]
    T_s = 45;
    T_av = (T_s + T_inf) / 2;
    k = 0.0265; %[W/mK]
    v = 1.62e-05; %[m^2/s]
    Pr = 0.711;
    beta = 1/(T_av + 273); %[1/K]
    
    %% Part a)
    
    % Guessing Value
    T_guess = 46.315;
    
    % Calculate Grashof Number
    Gr_ = Gr(g,beta,T_guess,T_inf,0.2,v);
    
    % Calculate Rayleigh Number
    Ra = Gr_ * Pr;
    
    % Calculate Nusselt Number
    Nu_ = Nu_VP(Ra,Pr);
    
    % Find convection heat coefficient
    h = Nu_ * k / 0.2; %[W/m^2K]
    
    % Heat Balance Equation (with radiation)
    %
    % Q_dis = Q_conv + Q_rad
    
    Q_trans_a = h*A_pcb*(T_guess - T_inf) + epsilon*sigma*A_pcb*((T_guess+273)^4 - (T_inf+273)^4);
    
    % After iterating across a couple of guesses to get the surface
    % temperature such that the heat leaving through convection and
    % radiation would be equal to the power dissapation, the resulting
    % surface temperature was:
    fprintf('Problem 2: \n \n')
    fprintf('Part a) \n')
    fprintf('Surface Temperature [deg C]: \n')
    disp(T_guess)
    
    
    %% Part b)
    
    % New Temperature Guess
    T_guess = 41.8; %[deg C]
    
    % Characteristic Length
    L_c = A_pcb / (2*0.15 + 2*0.2);
    
    % Grashof Number
    Gr_ = Gr(g,beta,T_guess,T_inf,L_c,v);
    
    % Rayleigh Number
    Ra_ = Gr_ * Pr;
    
    % Nusselt Number
    Nu_ = 0.59*Ra_^0.25;
    
    % Find convection heat coefficient
    h = Nu_ * k / L_c; %[W/m^2K]
    
    % Heat Balance Equation (with radiation)
    %
    % Q_dis = Q_conv + Q_rad
    
    Q_trans_b = h*A_pcb*(T_guess - T_inf) + epsilon*sigma*A_pcb*((T_guess+273)^4 - (T_inf+273)^4);
    
    % After iterating across a couple of guesses to get the surface
    % temperature such that the heat leaving through convection and
    % radiation would be equal to the power dissapation, the resulting
    % surface temperature was:
    fprintf('Part b) \n')
    fprintf('Surface Temperature [deg C]: \n')
    disp(T_guess)
    
    %% Part c)
    
    % New Temperature Guess
    T_guess = 50.1; %[deg C]
    
    % Characteristic Length
    L_c = A_pcb / (2*0.15 + 2*0.2);
    
    % Grashof Number
    Gr_ = Gr(g,beta,T_guess,T_inf,L_c,v);
    
    % Rayleigh Number
    Ra_ = Gr_ * Pr;
    
    % Nusselt Number
    Nu_ = 0.27*Ra_^0.25;
    
    % Find convection heat coefficient
    h = Nu_ * k / L_c; %[W/m^2K]
    
    % Heat Balance Equation (with radiation)
    %
    % Q_dis = Q_conv + Q_rad
    
    Q_trans_c = h*A_pcb*(T_guess - T_inf) + epsilon*sigma*A_pcb*((T_guess+273)^4 - (T_inf+273)^4);
    
    % After iterating across a couple of guesses to get the surface
    % temperature such that the heat leaving through convection and
    % radiation would be equal to the power dissapation, the resulting
    % surface temperature was:
    fprintf('Part c) \n')
    fprintf('Surface Temperature [deg C]: \n')
    disp(T_guess)


%% Problem 20-60
%
%

    %% Given
    eff = 0.1;
    d = 0.08; %[m]
    A = pi*d^2; %[m^2]
    P = 60; %[W]
    T_inf = 25; %[deg C]
    epsilon = 0.9;
    
    %% Solve
    
    % Make a temperature guess
    T_guess = 163; %[deg C]
    
    %% Properties of Surrounding Air (Table A-22)
    T_av = (T_guess + T_inf) / 2; %[deg C]
    k = 0.03095; %[W/mK]
    v = 2.306e-05; %[m^2/s]
    Pr = 0.7111;
    beta = 1/(T_av + 273);
    
    %% Calculate Convection & Radiation
    
    % Characteristic Length
    L_c = d;
    
    % Grashof Number
    Gr_ = Gr(g,beta,T_guess,T_inf,L_c,v);
    
    % Rayleigh Number
    Ra_ = Gr_ * Pr;
    
    % Nusselt Number
    Nu_ = 2 + (0.589*Ra_)^0.25/(1 + (0.469/Pr)^(9/16))^(4/9);
    
    % Find convection Heat coefficient
    h = Nu_ * k / L_c; %[W/m^2K]
    
    % Heat Balance Equation (with radiation)
    %
    % Q_dis = Q_conv + Q_rad
    
    Q_in = P * (1 - eff);
    
    Q_transfer = h*A*(T_guess - T_inf) + epsilon*sigma*A*((T_guess+273)^4 - (T_inf+273)^4);
    
    % After iterating across a couple of guesses to get the surface
    % temperature such that the heat leaving through convection and
    % radiation would be equal to the power absorption, the resulting
    % surface temperature was:
    fprintf('Problem 20-60: \n \n')
    fprintf('Surface Temperature [deg C]: \n')
    disp(T_guess)




