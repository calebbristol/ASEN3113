%% ASEN 3113 Lab 01
%
% Group Members:
%  -Caleb Bristol
%  -Jared Seefried
%  -Henri Wessels
%  -Sophia Trussel
%  -Justin Penderson
%
% Date: 1/12/22
%      


%% Housekeeping

clc
clear
close all;


%% Read In Data
data_8 = load('Data/Group18_Temp8');
data_10 = load('Data/Group18_Temp10');
data_12 = load('Data/Group18_Temp12');


%% Obtain RPM 
freq_8 = freqfinder(data_8);
freq_10 = freqfinder(data_10);
freq_12 = freqfinder(data_12);

freq_8_rpm = freq_8(end) * 60;
freq_10_rpm = freq_10(end) * 60;
freq_12_rpm = freq_12(end) * 60;


%% Obtain Average Temperature Difference
delt_8 = mean(data_8(:,7) - data_8(:,6));
delt_10 = mean(data_10(:,7) - data_10(:,6));
delt_12 = mean(data_12(:,7) - data_12(:,6));


%% Obtain Volume Measurements over Time


    %% Known Dimensions
    
    
    %% Read in Data
    [SH_ang_84] = readtable('Data/CAD/84_Small_Hole_Disp.csv');
    
    
    %% Find Frequencies 
    freq_84 = freqfinder([SH_ang_84.Var2 SH_ang_84.Var3]);
    
    
    %% Volume Calculations

    % bottom chamber
    d = 144;            % diameter, mm
    h = 21;             % height, mm
    b_cham_vol = h*pi*(d/2)^2;

    % bottom disk
    d_bd = 140;         % diameter, mm
    h_bd = 11;          % height, mm
    b_disk_vol = h_bd*pi*(h_bd/2)^2;

    % bottom volume
    b_vol = b_cham_vol - b_disk_vol;

    % small chamber max displacement
    h_sc_max = 25;
    d_sc = 14.92;
    t_vol_max = h_sc_max*pi*(d_sc/2)^2;

    % extra small chamber volume
    h_xtra_t = 6.5;
    xtra_vol = h_xtra_t*pi*(d_sc/2)^2;

    % max volume (mm^3)
    max_vol = b_vol + t_vol_max + xtra_vol;

    % min volume (mm^3)
    min_vol = b_vol + xtra_vol;

    % delta volume
    delta_vol = max_vol - min_vol;

    
    %% 


%% Functions 

function freq = freqfinder(data)
time = data(:,1);
p = data(:,2);

 % Calculate Frequency
    %
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    yidx = zci(p);
    for i = 2:length(yidx)-1
        tc(i) = interp1(p(yidx(i)+[-1:1]), time(yidx(i)+[-1:1]), 0, 'linear');  % ‘Exact’ Zero Crossings
    end
    frequency = 1./(2*gradient(tc));
    
    for i = 1:length(tc)
        if tc(i) > 120
            index_end = i;
            tc = tc(1:index_end);
            frequency = frequency(1:index_end);
            break;
        end
    end
    
    freq = frequency;
end