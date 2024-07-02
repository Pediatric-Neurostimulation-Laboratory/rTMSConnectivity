%% This script visualizes the change of wPLI connectivity after real vs.
% sham rTMS in all ROIs
%
% Author: Xiwei She, Ph.D.

clear; clc;

addpath('toolbox')

% Specify the results to check
dataPath = 'results\';

% Define subjects
subjList = {'example'}; % All subjects

s1 = {'example'}; % List all participants who had Right hemisphere rTMS
s2 = {}; % List all participants who had Left hemisphere rTMS

% Hyper-parameters
method = 'wpli_debiased';

timeRange = '10 to 1000';

freq = 'BETA';

%% This now loads all the data and creates a matrix of all data from all subjects, organized into ipsi-contra
for f = 1:length(subjList)

    % Load PRE rTMS and POST rTMS csv files
    sessionName = 'ses-2_pre_sham';
    PreSham = readmatrix([dataPath, (char(subjList(f))), '_', sessionName, '_', method, '_', freq, '_', timeRange, '.csv']);
    sessionName = 'ses-2_post_sham';
    PostSham = readmatrix([dataPath, (char(subjList(f))), '_', sessionName, '_', method, '_', freq, '_', timeRange, '.csv']);
    chanLocs_Sham = load([dataPath, (char(subjList(f))), '_', sessionName, 'channelLocations.mat']);

    sessionName = 'ses-1_pre_real';
    PreReal = readmatrix([dataPath, (char(subjList(f))), '_', sessionName, '_', method, '_', freq, '_', timeRange, '.csv']);
    sessionName = 'ses-1_post_real';
    PostReal = readmatrix([dataPath, (char(subjList(f))), '_', sessionName, '_', method, '_', freq, '_', timeRange, '.csv']);
    chanLocs_Real = load([dataPath, (char(subjList(f))), '_', sessionName, 'channelLocations.mat']);

    % Calculate Difference between Pre and Post Connectivity Matrices
    rTMS_Change_Real = PostReal - PreReal;

    rTMS_Change_Sham = PostSham - PreSham;

    rTMS_Real_v_Sham = rTMS_Change_Real - rTMS_Change_Sham;

    % Calculate Regional wPLI by averaging channel wPLI from 63*63 --> 9*9
    % Calculate X for MakeNodalMatrix: 1=left rTMS; 2=right rTMS;
    tf = strcmp(subjList{f}, s1);
    tf2 = strcmp(subjList{f}, s2);
    if mean(tf)>0
        x=2;
    elseif mean(tf2)>0
        x=1;
    end

    Nodal_Change_Real{f} = calculateRegionalConnectivity(rTMS_Change_Real, x, chanLocs_Real.chanLocs);
    Nodal_Change_Sham{f} = calculateRegionalConnectivity(rTMS_Change_Sham, x, chanLocs_Sham.chanLocs);
    Nodal_Change_Comparison{f} = calculateRegionalConnectivity(rTMS_Real_v_Sham, x, chanLocs_Real.chanLocs);

    Nodal_PreReal{f} = calculateRegionalConnectivity(PreReal, x, chanLocs_Real.chanLocs);
    Nodal_PostReal{f} = calculateRegionalConnectivity(PostReal, x, chanLocs_Real.chanLocs);
    Nodal_PreSham{f} = calculateRegionalConnectivity(PreSham, x, chanLocs_Sham.chanLocs);
    Nodal_PostSham{f} = calculateRegionalConnectivity(PostSham, x, chanLocs_Sham.chanLocs);


    % Reduce Data from 9*9 --> 6*6 Only consider Frontal Motor and Temporal
    temp = Nodal_Change_Real{f}; Nodal_Change_Real{f} = temp(1:6, 1:6);
    temp = Nodal_Change_Sham{f}; Nodal_Change_Sham{f} = temp(1:6, 1:6);
    temp = Nodal_Change_Comparison{f}; Nodal_Change_Comparison{f} = temp(1:6, 1:6);

    temp = Nodal_PreReal{f}; Nodal_PreReal{f} = temp(1:6, 1:6);
    temp = Nodal_PostReal{f}; Nodal_PostReal{f} = temp(1:6, 1:6);
    temp = Nodal_PreSham{f}; Nodal_PreSham{f} = temp(1:6, 1:6);
    temp = Nodal_PostSham{f}; Nodal_PostSham{f} = temp(1:6, 1:6);

end

ROIs = {'I Frontal', 'C Frontal', 'I Central', 'C Central', 'I Temporal', 'C Temporal'};
ROIs_vis = {'C-T', 'I-T', 'C-C', 'I-C', 'C-F', 'I-F'};

% Real rTMS
figure('Position', [1, 50, 1800, 600], 'Name', [subjList{f}, ' Real rTMS'])
colorLimit_redblue_1 = -0.15;
colorLimit_redblue_2 = 0.15;
t = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
nexttile(t)
hm_2 = heatmap(ROIs, ROIs, Nodal_PreReal{f}); hm_2.Colormap = gray; hm_2.CellLabelFormat = '%0.2g';
hm_2.ColorLimits = [0 0.15];
title('Connectivity Pre (Real)');
nexttile(t)
hm_3 = heatmap(Nodal_PostReal{f}); hm_3.Colormap = gray; hm_3.CellLabelFormat = '%0.2g';
hm_3.ColorLimits = [0 0.15];
title('Connectivity Post (Real)');
nexttile(t)
hm_6 = heatmap(Nodal_Change_Real{f});
hm_6.Colormap = redwhiteblue(colorLimit_redblue_1, colorLimit_redblue_2);
hm_6.ColorLimits = [colorLimit_redblue_1 colorLimit_redblue_2]; hm_6.CellLabelFormat = '%0.2g';
title('Connectivity Change Post - Pre');

% Sham rTMS
figure('Position', [1, 50, 1800, 600], 'Name', [subjList{f}, ' Sham rTMS'])
t = tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
nexttile(t)
hm_2 = heatmap(ROIs, ROIs, Nodal_PreSham{f}); hm_2.Colormap = gray; hm_2.CellLabelFormat = '%0.2g';
hm_2.ColorLimits = [0 0.15];
title('Connectivity Pre (Sham)');
nexttile(t)
hm_3 = heatmap(Nodal_PostSham{f}); hm_3.Colormap = gray; hm_3.CellLabelFormat = '%0.2g';
hm_3.ColorLimits = [0 0.15];
title('Connectivity Post (Sham)');
nexttile(t)
hm_6 = heatmap(Nodal_Change_Sham{f});
hm_6.Colormap = redwhiteblue(colorLimit_redblue_1, colorLimit_redblue_2);
hm_6.ColorLimits = [colorLimit_redblue_1 colorLimit_redblue_2]; hm_6.CellLabelFormat = '%0.2g';
title('Connectivity Change Post - Pre');
