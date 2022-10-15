%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ASEN 3113 - 18.19
%  Author: Caleb Bristol
%  Date: 4 April, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all;

%% Given
L = 0.02; %[m]
k = 21; %[W/mK]
rho = 8000; %[kg/m^3]
c_p = 570; %[J/kgK]
T_i = 18; %[deg C]
T_a = 950; %[deg C]
h = 150; %[W/m^2K]
d = 3; %[m]
v = 0.005:0.001:0.06; %[m/s]

%% Solve

    %% Verify Lumped Analysis
    L_c = L/2;
    
    Bi = h*L_c/k;
    
    fprintf('Biot Number: \n')
    disp(Bi)
    
    %% Run Simulation
    b = h/(rho*c_p*L_c);
    t_ = d./v;
    
    T = @(t) T_a + (T_i-T_a).*exp(-(b*t));
    
    T_ = T(t_);
    
    %% Plotting 
    figure()
    plot(v.*1000,T_,'LineWidth',2); hold on
    xlabel('Velocity [mm/s]')
    ylabel('Final Temperature [^oC]')
    title('Final Plate Temperature Against Plate Velocity')
    set(gca,'Fontsize',14)
    grid on; grid minor;
    hold off