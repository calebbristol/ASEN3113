%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ASEN 3113: Design Lab
%  Group Members:
%  -
%
%  Date:
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all;


%% Constants
R_e = 6378; %[km]
R_e_o = 1 * 1.495978707e08; %[km]
R_s = R_e_o / 215; %[km]
a = R_e_o;
R_sc_o = 1 * 42164; %[km]
P_sc = 60*60*24; %[s]
inc = 22.5; %[deg]
e = 0.01671; %[]
epsilon = 1;
sigma = 5.67e-08; %[W/m^2K^4]
T_sun = 5778; %[K]
mu_s = 1.32712440018e11; %[km^3/s^2]
mu_e = 3.986004418e05; %[km^3/s^2]


%% Solar Radiation
%
% This section calculates the solar radiation energy available to the
% spacecraft. First, Earth's distance to the sun is calculated, then the 
G_s = epsilon*sigma*(T_sun)^4; %[W/m^2]
E_solar = @(r) R_s^2./r.^2.*G_s; %[W/m^2]


    %% Distance From Sun
    % Perihelion occurs on January 4
    
    % Orbital Mechanics 
    p = a*(1-e^2);
    n = sqrt(mu_s/a^3);
    P = 2*pi/n;
    M = @(t) 2*pi/P*t;
    t = 0:60*60*24:P;
    M_ = M(t);
    
    % Preallocate Vectors
    E_ = zeros(length(t),1);
    f_ = zeros(length(t),1);
    r_ = zeros(length(t),1);
    
    % Run Kepler's time of flight integration
    for i = 1:length(t)
        [E_(i),f_(i)] = keptof(2*pi*i/length(t),e,M_(i),1e-09);
    end
    
    % Calculate Radius at Every Day of Year
    r = a.*(1-e.*cos(E_));
    
    % Rotate for January 4 offset
    r = circshift(r,3);
    

    %% Summer Solstice
    % Summer Solstice occurs on June 21
    
    % Calculate Earth Radius
    r_ss = r(31+28+31+30+31+21);
    
    % Add effects of satellite orbit
    r_sc_ss = r_ss + R_sc_o*-cos(linspace(0,2*pi,P_sc)).*cosd(inc);
    
    % Calcualte Energy
    E_ss = E_solar(r_sc_ss);

    
    %% Winter Solstice
    % Winter Solstice occurs on January 21
    
    % Calculate Earth Radius
    r_ws = r(21);
    
    % Add effects of satellite orbit
    r_sc_ws = r_ws + R_sc_o*-cos(linspace(0,2*pi,P_sc)).*cosd(inc);
    
    % Calculate Energy
    E_ws = E_solar(r_sc_ws);
    
    
    %% Equinox
    % Equinox occurs on March 20 or September 22
    r_eq = r(31+28+20);
    
    % Add effects of satellite orbit
    r_sc_eq = r_eq + R_sc_o*-cos(linspace(0,2*pi,P_sc));
    
    % Calculate Energy
    E_eq = E_solar(r_sc_eq);
    
    % Check for Eclipse
    eclipse_mult = ones(1,length(E_eq)); % This is used for IR backload
    
    for i = 1:P_sc
        if abs(R_sc_o*sin(2*pi*i/P_sc)) < R_e && -cos(2*pi*i/P_sc) > 0
            E_eq(i) = 0;
            eclipse_mult(i) = 11/75.5;
        end
    end

    
%% Incident Radiation
%
%


    %% Summer Solstice
    Rad_ss = max(E_ss.*sin(linspace(0,2*pi,P_sc)).*cosd(inc),0);
    
    
    %% Winter Solstice
    Rad_ws = max(E_ws.*sin(linspace(0,2*pi,P_sc)).*cosd(inc),0);
    
    
    %% Equinox
    Rad_eq = max(E_eq.*sin(linspace(0,2*pi,P_sc)),0);
  
    
%% Plotting
%
%


    %% Summer Solstice
    figure()
    plot(1:P_sc,Rad_ss,'LineWidth',2); hold on
    
    %% Winter Solstice
    plot(1:P_sc,Rad_ws,'LineWidth',2);
    
    %% Equinox
    plot(1:P_sc,Rad_eq,'LineWidth',2);
    
    title('Sun Irradiation')
    xlabel('Time [s]')
    ylabel('Irradiation G [W/m^2]')
    legend('Summer Solstice','Winter Solstice','Equinox')
    grid on; grid minor;
    hold off
    
    
%% Radiator Area
%
%


    %% Energy Balance Equation
    %
    % Q_sun + Q_heater + Q_instr + Q_IR = Q_rad
    % Q_sun + Q_IR - Q_rad = -Q_heater - Q_instr
    %
    % E_sun(A) + \sigma + Heater + 20/0 W + IR(A) = 0
    alpha_sol = 0.2;
    e_rad = 0.85;
    alpha_rad = 0.85;
    
    A = (-0 - 20) / (alpha_sol*max(E_eq) + alpha_rad*75.5 - e_rad*sigma*(303^4-2.7^4));
    
    A_1 = A;
    
    
%% Total Radiation Flux
%
%


    %% Function a la Irradiance & Area
    Flux = @(Area,G_s,alpha) Area .* G_s .* alpha;
    
    
    %% Summer Solstice
    figure()
    plot(1:P_sc,Flux(A,Rad_ss,alpha_sol) + alpha_rad*63*A - e_rad*sigma*A*303^4,'LineWidth',2); hold on
    
    %% Winter Solstice
    plot(1:P_sc,Flux(A,Rad_ws,alpha_sol) + alpha_rad*88*A - e_rad*sigma*A*303^4,'LineWidth',2);
    
    %% Equinox
    plot(1:P_sc,Flux(A,Rad_eq,alpha_sol)+ alpha_rad*75.5*A*eclipse_mult - e_rad*sigma*A*303^4,'LineWidth',2);
    
    title('Total Radiation Flux')
    xlabel('Time [s]')
    ylabel('Radiation Flux [W]')
    legend('Summer Solstice','Winter Solstice','Equinox')
    grid on; grid minor;
    hold off
    
    
%% Total Solar Radiation Flux
%
%

    %% Summer Solstice
    figure()
    plot(1:P_sc,Flux(1,Rad_ss,alpha_sol),'LineWidth',2); hold on
    
    %% Winter Solstice
    plot(1:P_sc,Flux(1,Rad_ws,alpha_sol),'LineWidth',2);
    
    %% Equinox
    plot(1:P_sc,Flux(1,Rad_eq,alpha_sol),'LineWidth',2);
    
    title('Total Solar Radiation Flux')
    xlabel('Time [s]')
    ylabel('Radiation Flux [W/m^2]')
    legend('Summer Solstice','Winter Solstice','Equinox')
    grid on; grid minor;
    hold off
    
    
    
%% Operating Temperature
%
% 

    
    %% Summer Solstice
    T_ss_op = ((Flux(A,Rad_ss,alpha_sol) + alpha_rad*63*A + 20)./(e_rad.*sigma.*A)).^0.25;
    T_ss_surv = ((Flux(A,Rad_ss,alpha_sol) + alpha_rad*63*A)./(e_rad.*sigma.*A)).^0.25;
    
    %% Winter Solstice
    T_ws_op = ((Flux(A,Rad_ws,alpha_sol) + alpha_rad*88*A + 20)./(e_rad.*sigma.*A)).^0.25;
    T_ws_surv = ((Flux(A,Rad_ws,alpha_sol) + alpha_rad*88*A)./(e_rad.*sigma.*A)).^0.25;
    
    %% Equinox
    T_eq_op = ((Flux(A,Rad_eq,alpha_sol) + alpha_rad*75.5*A*eclipse_mult + 20)./(e_rad.*sigma.*A)).^0.25;
    T_eq_surv = ((Flux(A,Rad_eq,alpha_sol) + alpha_rad*75.5*A*eclipse_mult)./(e_rad.*sigma.*A)).^0.25;
    
    
%% Heater Power
%
%


    %% Summer Solstice
    P_ss_op = max(-Flux(A,Rad_ss,alpha_sol) - alpha_rad*63*A + e_rad*sigma*A*303^4 - 20,0);
    P_ss_surv = max(-Flux(A,Rad_ss,alpha_sol) - alpha_rad*63*A + e_rad*sigma*A*233^4,0);
    
    %% Winter Solstice
    P_ws_op = max(-Flux(A,Rad_ws,alpha_sol) - alpha_rad*88*A + e_rad*sigma*A*303^4 - 20,0);
    P_ws_surv = max(-Flux(A,Rad_ws,alpha_sol) - alpha_rad*88*A + e_rad*sigma*A*233^4,0);
    
    %% Equinox
    P_eq_op = max(-Flux(A,Rad_eq,alpha_sol) - alpha_rad*75.5*A*eclipse_mult + e_rad*sigma*A*303^4 - 20,0);
    P_eq_surv = max(-Flux(A,Rad_eq,alpha_sol) - alpha_rad*75.5*A*eclipse_mult + e_rad*sigma*A*233^4,0);
    
    
%% Plotting
%
%


    %% Summer Solstice
    figure()
    yyaxis left
    plot(1:P_sc,T_ss_op - 273); hold on
    plot(1:P_sc,T_ss_surv - 273)
    ylabel('Temperature [^oC]')
    yyaxis right
    plot(1:P_sc,P_ss_op)
    plot(1:P_sc,P_ss_surv)
    title('Summer Solstice')
    xlabel('Time [s]')
    ylabel('Power [W]')
    ylim([0 160])
    legend('Operational Unheated Temperature','Survival Unheated Temperature','Operational Heater Power','Survival Heater Power')
    hold off
    
    %% Winter Solstice
    figure()
    yyaxis left
    plot(1:P_sc,T_ws_op - 273); hold on
    plot(1:P_sc,T_ws_surv - 273)
    ylabel('Temperature [^oC]')
    yyaxis right
    plot(1:P_sc,P_ws_op)
    plot(1:P_sc,P_ws_surv)
    title('Winter Solstice')
    xlabel('Time [s]')
    ylabel('Power [W]')
    ylim([0 160])
    legend('Operational Unheated Temperature','Survival Unheated Temperature','Operational Heater Power','Survival Heater Power')
    hold off
    
    %% Equinox
    figure()
    yyaxis left
    plot(1:P_sc,T_eq_op - 273); hold on
    plot(1:P_sc,T_eq_surv - 273)
    ylabel('Temperature [^oC]')
    yyaxis right
    plot(1:P_sc,P_eq_op)
    plot(1:P_sc,P_eq_surv)
    title('Equinox')
    xlabel('Time [s]')
    ylabel('Power [W]')
    ylim([0 160])
    legend('Operational Unheated Temperature','Survival Unheated Temperature','Operational Heater Power','Survival Heater Power')
    hold off
   
    
    
%% New Coating
%
% Case 1: Aluminum Silicate/ Potassium Silicate

    %% Radiation Properties
    alpha_sol = 0.14;
    e_rad = 0.94;
    alpha_rad = 0.94;
    
    A = (-0 - 20) / (alpha_sol*max(E_eq) + alpha_rad*75.5 - e_rad*sigma*(303^4 - 2.7^4));
    
    A_2 = A;
    
    
    %% Radiation Flux
    figure()
    plot(1:P_sc,Flux(A,Rad_eq,alpha_sol)+ alpha_rad*75.5*A*eclipse_mult - e_rad*sigma*A*303^4,'LineWidth',2);
    
    title('Total Radiation Flux, Improved Coating')
    xlabel('Time [s]')
    ylabel('Radiation Flux [W]')
    legend('Equinox')
    grid on; grid minor;
    hold off
    
    
    %% Solar Radiation Flux
    figure()
    plot(1:P_sc,Flux(1,Rad_eq,alpha_sol),'LineWidth',2);
    
    title('Total Solar Radiation Flux, Improved Coating')
    xlabel('Time [s]')
    ylabel('Radiation Flux [W/m^2]')
    legend('Equinox')
    grid on; grid minor;
    hold off
    
    %% Temperature and Power
    T_eq_op = ((Flux(A,Rad_eq,alpha_sol) + alpha_rad*75.5*A*eclipse_mult + 20)./(e_rad.*sigma.*A)).^0.25;
    T_eq_surv = ((Flux(A,Rad_eq,alpha_sol) + alpha_rad*75.5*A*eclipse_mult)./(e_rad.*sigma.*A)).^0.25;
    
    P_eq_op = max(-Flux(A,Rad_eq,alpha_sol) - alpha_rad*75.5*A*eclipse_mult + e_rad*sigma*A*303^4 - 20,0);
    P_eq_surv = max(-Flux(A,Rad_eq,alpha_sol) - alpha_rad*75.5*A*eclipse_mult + e_rad*sigma*A*233^4,0);

    figure()
    yyaxis left
    plot(1:P_sc,T_eq_op - 273); hold on
    plot(1:P_sc,T_eq_surv - 273)
    ylabel('Temperature [^oC]')
    yyaxis right
    plot(1:P_sc,P_eq_op)
    plot(1:P_sc,P_eq_surv)
    title('Equinox, Improved Coating')
    xlabel('Time [s]')
    ylabel('Power [W]')
    ylim([0 160])
    legend('Operational Unheated Temperature','Survival Unheated Temperature','Operational Heater Power','Survival Heater Power')
    hold off
    
    
%% New Coating
%
% Case 2: 

    %% Radiation Properties
    alpha_sol = 0.14;
    e_rad = 0.94;
    alpha_rad = 0.94;
    
    A = (-0 - 20) / (alpha_sol*max(E_eq) + alpha_rad*75.5 - e_rad*sigma*(303^4 - 2.7^4));
    
    A_3 = A;
    
    
    %% Radiation Flux
    figure()
    plot(1:P_sc,Flux(A,Rad_eq,alpha_sol)+ alpha_rad*75.5*A*eclipse_mult - e_rad*sigma*A*303^4,'LineWidth',2);
    
    title('Total Radiation Flux, Improved Coating')
    xlabel('Time [s]')
    ylabel('Radiation Flux [W]')
    legend('Equinox')
    grid on; grid minor;
    hold off
    
    
    %% Solar Radiation Flux
    figure()
    plot(1:P_sc,Flux(1,Rad_eq,alpha_sol),'LineWidth',2);
    
    title('Total Solar Radiation Flux, Improved Coating')
    xlabel('Time [s]')
    ylabel('Radiation Flux [W/m^2]')
    legend('Equinox')
    grid on; grid minor;
    hold off
    
    
    %% Temperature and Power
    T_eq_op = ((Flux(A,Rad_eq,alpha_sol) + alpha_rad*75.5*A*eclipse_mult + 20)./(e_rad.*sigma.*A)).^0.25;
    T_eq_surv = ((Flux(A,Rad_eq,alpha_sol) + alpha_rad*75.5*A*eclipse_mult)./(e_rad.*sigma.*A)).^0.25;
    
    P_eq_op = max(-Flux(A,Rad_eq,alpha_sol) - alpha_rad*75.5*A*eclipse_mult + e_rad*sigma*A*303^4 - 20,0);
    P_eq_surv = max(-Flux(A,Rad_eq,alpha_sol) - alpha_rad*75.5*A*eclipse_mult + e_rad*sigma*A*233^4,0);

    figure()
    yyaxis left
    plot(1:P_sc,T_eq_op - 273); hold on
    plot(1:P_sc,T_eq_surv - 273)
    ylabel('Temperature [^oC]')
    yyaxis right
    plot(1:P_sc,P_eq_op)
    plot(1:P_sc,P_eq_surv)
    title('Equinox, Improved Coating')
    xlabel('Time [s]')
    ylabel('Power [W]')
    ylim([0 160])
    legend('Operational Unheated Temperature','Survival Unheated Temperature','Operational Heater Power','Survival Heater Power')
    hold off
    

%% Functions

function [E,f] = keptof(E_0,e,M,tol)
f_calc = @(E) 2*atan(sqrt((1+e)/(1-e)) * tan(E/2));
ratio = 1;
E_ = E_0;
%f_ = f_calc(E_0);

%i = 2;
    while ratio > tol
        ratio = (E_ - e*sin(E_) - M)/(1 - e*cos(E_));
        E_ = E_ - ratio;
        %f_(i) = f_calc(E_(i));
        %i = i+1;
    end
%iteration = 1:i-1;
E = E_;
f = f_calc(E_);
end