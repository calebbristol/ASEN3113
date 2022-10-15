%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ASEN 3113: Homework 8
%  Author: Caleb Bristol
%
%  Date: 13 April, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all;

%% Chapter 21: Problem 21


    %% Given
    T_sun = 5870; %[K]
    lambda = 0.01:0.01:1000; %[um]
    
    %% Emissive Power Function
    %
    % Utilizing Planck's Law
    C_1 = 3.74177e08; %[W um^4/m^2]
    C_2 = 1.43878e04; %[um K]
    
    E_b_lambda = @(lambda,T) C_1./(lambda.^5.*(exp(C_2./(lambda.*T))-1));
    
    %% Plotting
    figure()
    plot(lambda,E_b_lambda(lambda,T_sun),'r','LineWidth',2); hold on
    title('Spectral Blackbody Emissive Power vs. Wavelength (Sun)')
    xlabel('\lambda [um]')
    ylabel('E_b_\lambda [W um/m^2]')
    xlim([0.01 1000])
    ylim([10e-06 10e08])
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    grid on; grid minor;
    hold off
    
    figure()
    plot(lambda,E_b_lambda(lambda,T_sun),'r','LineWidth',2); hold on
    title('Spectral Blackbody Emissive Power vs. Wavelength (Sun)')
    xlabel('\lambda [um]')
    ylabel('E_b_\lambda [W um/m^2]')
    xlim([0.01 1000])
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    grid on; grid minor;
    hold off
    
    %% Discussion
    %
    % The plot is almost entirely linear when viewed with linear axis
    % scaling. When put in logarithmic scaling, it appears as a much
    % smoother curve, though much more level than the respective figure
    % shown in the textbook. The textbook cuts off the bottom end of the
    % scale significantly, so it hides the fact that there is almost no
    % short wavelength waves leaving the sun (though that is implied), and
    % there is still a substantial amount of waves coming from the infrared
    % regions. However, there is a significant amount of radiation coming
    % through the visible spectrum, and it is exponentially higher than
    % that around the rest of the spectrum.