% Main script for QRS detection using morphological operations

% Clear workspace and command window
close all 
clear all
clc

% Set the main directory 
mainDir = fileparts(which('main_morphological.m'));
if isempty(mainDir)
    error('main_morphological.m not found in the MATLAB path. Please ensure it exists and is in the MATLAB path.');
end

% Set the data path
dataPath = fullfile(mainDir, 'mit-bih-arrhythmia-database-1.0.0');

% Ensure the data directory exists
if ~exist(dataPath, 'dir')
    error('Data directory not found. Please ensure the MIT-BIH Arrhythmia Database is in the correct location.');
end

% Add the WFDB Toolbox to the MATLAB path
addpath(genpath(mainDir));

% Set the current directory to the data path
cd(dataPath);

% Specify the record name (e.g., 108, 109, 119, 212)
record_name = '108';  

% Read ECG signal
try
    [signal, Fs, tm] = rdsamp(record_name, 1);
    disp(['Successfully read ' num2str(length(signal)) ' samples']);
catch ME
    error('Error reading record: %s. Error message: %s', record_name, ME.message);
end

% Read annotations
try
    [ann, anntype, ~, ~, ~] = rdann(record_name, 'atr');
    disp(['Successfully read ' num2str(length(ann)) ' annotations']);
catch ME
    warning('Could not read annotations: %s', ME.message);
    ann = [];
end

% Detect QRS complexes using the morphological approach
[qrs_peaks, sensitivity, positive_predictivity] = detect_qrs_morphological(signal, Fs, ann);

% Create time vector
t = (0:length(signal)-1) / Fs;

% Plot ECG signal with detected QRS complexes and annotations
figure;
plot(t, signal);
hold on;

% Plot detected QRS complexes
plot(t(qrs_peaks), signal(qrs_peaks), 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% Plot annotations if available
if ~isempty(ann)
    plot(ann/Fs, signal(ann), 'g^', 'MarkerSize', 8, 'LineWidth', 2);
end

% Set plot properties
xlabel('Time (s)');
ylabel('Amplitude');
title(['ECG Signal with Detected QRS Complexes and Annotations - Record ', record_name]);
legend('ECG Signal', 'Detected QRS', 'Annotations');
grid on;

% Display performance metrics if available
if ~isnan(sensitivity) && ~isnan(positive_predictivity)
    fprintf('QRS Detection Performance:\n');
    fprintf('Sensitivity: %.2f%%\n', sensitivity * 100);
    fprintf('Positive Predictive Value: %.2f%%\n', positive_predictivity * 100);
else
    fprintf('Performance metrics not available.\n');
end

% Change the current directory back to the main directory
cd(mainDir);