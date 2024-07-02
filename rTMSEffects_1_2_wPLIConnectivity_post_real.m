%% This script is used for calculating the connectivity measured by wPLI
% from scalp EEG recorded from children with SeLECTS during a TMS session
%
% Author: Xiwei She, Ph.D.

clear;clc

%% Add necessary data & functions to the path
addpath('toolbox');

% Dowload and extract the eeglab and fieldtrip into the toolbox folder
% Replace the name/version accordingly
addpath('toolbox/eeglab2022.0');
addpath('toolbox/fieldtrip-20220104');
addpath('toolbox/fieldtrip-20220104/utilities');
eeglab();

%% Calculate wPLI Connectivity on Session 1 (Real rTMS in this example)

% Specify path
dataPath = 'exampleData';
inputFileName = 'example_ses-1_post_real.set';
outputFilePrefix = 'example_ses-1_post_real';

% Specify Selected Time Period
timeSpanList = {[1001, 2001]}; % 10 ~ 1000 ms
offsetList = {10};
timeSpanNameList = {'10 to 1000'};

% You can check other time-window by un-commenting below
% timeSpanList = {[1, 1001], [1001, 1512], [1001, 2002], [1001, 2500], [1001, 2002], [1001, 2500]};
% offsetList = {-10, 10, 10, 10, 490, 990};
% timeSpanNameList = {'-1000 to -10', '10 to 510', '10 to 1000', '10 to 1500', '500 to 1000', '1000 to 1500'};

for t = 1:length(timeSpanList)

    selectedTimeSpan = timeSpanList{t};
    timeSpanName = timeSpanNameList{t};
    offset = offsetList{t};

    %% Load file
    EEG = pop_loadset([dataPath, '/', inputFileName]);

    % Define and Select Time Range
    EEG = pop_select(EEG, 'point', selectedTimeSpan);

    %% Run Alto wPLI Connectivity
    ALTO_METHODS = {'wpli_debiased'};
    ALTO_FREQ = {'ALPHA','BETA','THETA'};
    ALTO_KERNEL = eye(length(EEG.chanlocs)); %identity matrix for Alto call

    out = ALTO_spTMSEEGConnectivity(EEG, ALTO_KERNEL, 'channel', ALTO_METHODS, ALTO_FREQ, offset);

    %% Save Results
    RESULTS_DIR = 'results';
    for f = 1:length(ALTO_FREQ)
        for m = 1:length(ALTO_METHODS)

            freq = ALTO_FREQ{f}; METHOD = ALTO_METHODS{m};
            writematrix(out.(freq).(METHOD),RESULTS_DIR+"\" + outputFilePrefix+"_" + METHOD+"_" + freq+"_" + timeSpanName + ".csv");
        end

    end

    if t ~= length(timeSpanList)
        clearvars EEG
    end
end
%% Save channel locations
chanLocs = EEG.chanlocs;
save(RESULTS_DIR+"\" + outputFilePrefix+"channelLocations.mat", "chanLocs")