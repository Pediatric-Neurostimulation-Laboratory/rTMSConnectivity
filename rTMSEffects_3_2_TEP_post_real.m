%% This script is used for calculating the TMS evoked potential (TEP) from 
% scalp EEG recorded from children with SeLECTS during a TMS session
%
% Author: Xiwei She

clear;clc

%% Add necessary data & functions to the path
addpath('toolbox');

% Dowload and extract the eeglab into the toolbox folder
% Replace the name/version accordingly
addpath('toolbox/eeglab2022.0');
eeglab();

%% Information that require user input

% Specify path
folderName = 'exampleData';
inputFileName = 'example_ses-1_post_real.set';
outputFilePrefix = 'example_ses-1_post_real';

% Specify the hemisphere of stimulation
spTMSHemisphere = 'right';

%% Load file
EEG = pop_loadset([folderName, '/', inputFileName]);

%% Hyper-parameters

% Calculate area under curve within the time window
timeMin = 10;
% timeMax = 80;
timeMax = 1000;

%% Local TEPs
if strcmp(spTMSHemisphere, 'left')
    ROI_electrodes.Ipsi_Frontal = {'Fp1', 'F3', 'F7', 'AF7', 'AF3', 'F1', 'F5'};
    ROI_electrodes.Contra_Frontal = {'F4', 'F8', 'Fp2', 'F6', 'AF8', 'AF4', 'F2'};

    ROI_electrodes.Ipsi_Motor = {'FC5', 'FC1', 'C3', 'CP5', 'CP1', 'FC3', 'C1', 'C5', 'CP3'};
    ROI_electrodes.Contra_Motor = {'CP6', 'CP2', 'C4', 'FC6', 'FC2', 'CP4', 'C6', 'C2', 'FC4'};

    ROI_electrodes.Ipsi_Temporal = {'FT9', 'T7', 'TP9', 'P7', 'FT7', 'TP7'};
    ROI_electrodes.Contra_Temporal = {'P8', 'TP10', 'T8', 'FT10', 'TP8', 'FT8'};

else
    ROI_electrodes.Ipsi_Frontal = {'F4', 'F8', 'Fp2', 'F6', 'AF8', 'AF4', 'F2'};
    ROI_electrodes.Contra_Frontal = {'Fp1', 'F3', 'F7', 'AF7', 'AF3', 'F1', 'F5'};

    ROI_electrodes.Ipsi_Motor = {'CP6', 'CP2', 'C4', 'FC6', 'FC2', 'CP4', 'C6', 'C2', 'FC4'};
    ROI_electrodes.Contra_Motor = {'FC5', 'FC1', 'C3', 'CP5', 'CP1', 'FC3', 'C1', 'C5', 'CP3'};

    ROI_electrodes.Ipsi_Temporal = {'P8', 'TP10', 'T8', 'FT10', 'TP8', 'FT8'};
    ROI_electrodes.Contra_Temporal = {'FT9', 'T7', 'TP9', 'P7', 'FT7', 'TP7'};
end

EEG = pop_tesa_tepextract(EEG, 'ROI', 'tepName', 'IpsiFrontal', 'elecs', ROI_electrodes.Ipsi_Frontal);
EEG = pop_tesa_tepextract(EEG, 'ROI', 'tepName', 'ContraFrontal', 'elecs', ROI_electrodes.Contra_Frontal);
EEG = pop_tesa_tepextract(EEG, 'ROI', 'tepName', 'IpsiMotor', 'elecs', ROI_electrodes.Ipsi_Motor);
EEG = pop_tesa_tepextract(EEG, 'ROI', 'tepName', 'ContraMotor', 'elecs', ROI_electrodes.Contra_Motor);
EEG = pop_tesa_tepextract(EEG, 'ROI', 'tepName', 'IpsiTemporal', 'elecs', ROI_electrodes.Ipsi_Temporal);
EEG = pop_tesa_tepextract(EEG, 'ROI', 'tepName', 'ContraTemporal', 'elecs', ROI_electrodes.Contra_Temporal);

TEP_y_Ipsi_Frontal = EEG.ROI.IpsiFrontal.tseries;
TEP_y_Contra_Frontal = EEG.ROI.ContraFrontal.tseries;
TEP_y_Ipsi_Motor = EEG.ROI.IpsiMotor.tseries;
TEP_y_Contra_Motor = EEG.ROI.ContraMotor.tseries;
TEP_y_Ipsi_Temporal = EEG.ROI.IpsiTemporal.tseries;
TEP_y_Contra_Temporal = EEG.ROI.ContraTemporal.tseries;

% Only consider TEPs within the time window
indices = timeMin <= EEG.times & EEG.times <= timeMax;

x_subset = EEG.times(indices);

% trapezoidal integration for AUC
aucTEP.Ipsi_Frontal = trapz(x_subset, abs(TEP_y_Ipsi_Frontal(indices))); 
aucTEP.Contra_Frontal = trapz(x_subset, abs(TEP_y_Contra_Frontal(indices))); 
aucTEP.Ipsi_Motor = trapz(x_subset, abs(TEP_y_Ipsi_Motor(indices))); 
aucTEP.Contra_Motor = trapz(x_subset, abs(TEP_y_Contra_Motor(indices))); 
aucTEP.Ipsi_Temporal = trapz(x_subset, abs(TEP_y_Ipsi_Temporal(indices))); 
aucTEP.Contra_Temporal = trapz(x_subset, abs(TEP_y_Contra_Temporal(indices))); 

% Peak analyses on local tep (ipsi motor)
peakLatency_pos = 60; % specify peak of interest
peakRange_pos = [40, 80];  % time range where looking for peak

peakLatency_neg = 100;
peakRange_neg = [80, 120];

EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'positive', peakLatency_pos, peakRange_pos, 'method' ,'largest', 'samples', 5, 'tepName', 'IpsiMotor');
EEG = pop_tesa_peakanalysis( EEG, 'ROI', 'negative', peakLatency_neg, peakRange_neg, 'method' ,'largest', 'samples', 5, 'tepName', 'IpsiMotor');


%% Visualize Results
tablePlot = 'on'; % on off

% Peak analysis on the stimulated site (ipsi motor)
peakAnalysesResults = pop_tesa_peakoutput(EEG, 'tepName', 'IpsiMotor', 'calcType', 'amplitude', 'winType', 'individual', 'averageWin', 5, 'fixedPeak', [], 'tablePlot', tablePlot);
pop_tesa_plot(EEG, 'tepType', 'ROI', 'tepName', 'IpsiMotor', 'CI', 'on'); set(gcf, 'Position', [50 558 560 420])

disp(' ')
disp(['TEP AUC (Ipsi Frontal): ', mat2str(aucTEP.Ipsi_Frontal)]);
disp(['TEP AUC (Contra Frontal): ', mat2str(aucTEP.Contra_Frontal)]);
disp(['TEP AUC (Ipsi Motor): ', mat2str(aucTEP.Ipsi_Motor)]);
disp(['TEP AUC (Contra Motor): ', mat2str(aucTEP.Contra_Motor)]);
disp(['TEP AUC (Ipsi Temporal): ', mat2str(aucTEP.Ipsi_Temporal)]);
disp(['TEP AUC (Contra Temporal): ', mat2str(aucTEP.Contra_Temporal)]);
disp(' ')

%% Save Results
RESULTS_DIR = 'results';
save(RESULTS_DIR + "\" + outputFilePrefix + "_TEPPeaks_" +  mat2str(timeMin) + "-" + mat2str(timeMax) + "ms.mat", "peakAnalysesResults", "aucTEP")