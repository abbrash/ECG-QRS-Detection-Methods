function [R_peaks] = detect_qrs_emd(ecg_signal, fs)
% QRS Complex Detection Algorithm using EMD and Non-linear Transformation
%
% Input:
%   ecg_signal: raw ECG signal
%   fs: sampling frequency in Hz
%
% Output:
%   R_peaks: indices of detected R peaks

% Define chunk size (e.g., 5 seconds of data)
chunk_size = 5 * fs;
signal_length = length(ecg_signal);
num_chunks = ceil(signal_length / chunk_size);

% Initialize output
R_peaks = [];

% Process signal in chunks to handle long recordings
for i = 1:num_chunks
    start_idx = (i-1)*chunk_size + 1;
    end_idx = min(i*chunk_size, signal_length);
    chunk = ecg_signal(start_idx:end_idx);

    % Step 1: Preprocessing - High-pass filter to remove baseline wander
    [b, a] = my_butter(4, 0.5/(fs/2), 'high');
    ecg_filtered = filter(b, a, chunk);

    % Step 2: Simplified EMD (extract single IMF)
    imf = extract_imf(ecg_filtered, fs);

    % Step 3 & 4: Non-linear Transformation and Moving Window Integration
    w = round(0.1 * fs); % Window size (100 ms)
    window = ones(1, w) / w;
    snt = non_linear_transform(imf);

    % Debugging output
    disp(['Size of snt: ', num2str(size(snt))]);
    disp(['Size of window: ', num2str(size(window))]);

    % Ensure snt is a row vector
    snt = snt(:)';
    si = conv(snt, window, 'same');

    % Step 5 & 6: Low-pass Filtering
    [b, a] = my_butter(1, 2/(fs/2), 'low');
    smooth_signal = filter(b, a, si);

    % Step 7: Peak Detection
    [~, chunk_peaks] = simple_findpeaks(smooth_signal, round(0.2*fs));

    % Adjust peak indices to global signal
    chunk_peaks = chunk_peaks + start_idx - 1;

    % Append to overall R_peaks
    R_peaks = [R_peaks, chunk_peaks];
end
end

function imf = extract_imf(signal, fs)
% Extract the first Intrinsic Mode Function (IMF) using simplified EMD
max_iterations = 10;
imf = signal;
for iter = 1:max_iterations
    [upper_env, lower_env] = envelope(imf, fs);
    mean_env = (upper_env + lower_env) / 2;
    prev_imf = imf;
    imf = imf - mean_env;
    if is_imf(imf)
        break;
    end
end
end

function [upper_env, lower_env] = envelope(signal, fs)
% Compute upper and lower envelopes of the signal
[~, max_idx] = simple_findpeaks(signal, round(0.2*fs));
[~, min_idx] = simple_findpeaks(-signal, round(0.2*fs));
if isempty(max_idx) || isempty(min_idx)
    upper_env = signal;
    lower_env = signal;
else
    upper_env = interp1(max_idx, signal(max_idx), 1:length(signal), 'spline', 'extrap');
    lower_env = interp1(min_idx, signal(min_idx), 1:length(signal), 'spline', 'extrap');
end
end

function result = is_imf(signal)
% Check if the signal satisfies IMF conditions
[~, max_idx] = simple_findpeaks(signal, 1);
[~, min_idx] = simple_findpeaks(-signal, 1);
result = abs(length(max_idx) + length(min_idx) - length(signal)) <= 1;
end

function [pks, locs] = simple_findpeaks(signal, minPeakDistance)
% Simple peak detection algorithm
pks = [];
locs = [];
signalLength = length(signal);
for i = 2 : (signalLength - 1)
    if (signal(i) > signal(i-1)) && (signal(i) > signal(i+1))
        if ~isempty(locs) && (i - locs(end) < minPeakDistance)
            if signal(i) > pks(end)
                pks(end) = signal(i);
                locs(end) = i;
            end
        else
            pks(end+1) = signal(i);
            locs(end+1) = i;
        end
    end
end
end

function snt = non_linear_transform(imf)
% Apply non-linear transformation to enhance QRS complexes
snt = zeros(size(imf));
for j = 3:length(imf)
    if (imf(j) * imf(j-1) > 0) && (imf(j) * imf(j-2) > 0)
        snt(j) = abs(imf(j) * imf(j-1) * imf(j-2));
    end
end
snt = snt(:)'; % Ensure output is a row vector
end