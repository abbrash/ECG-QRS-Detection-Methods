function [qrs_peaks, sensitivity, positive_predictivity] = detect_qrs_morphological(signal, fs, annotations)
    % QRS Detection using Morphological Operations
    %
    % This function detects QRS complexes in an ECG signal using morphological operations.
    %
    % Inputs:
    %   signal: Input ECG signal (multi-channel)
    %   fs: Sampling frequency of the ECG signal
    %   annotations: Ground truth annotations for performance evaluation
    %
    % Outputs:
    %   qrs_peaks: Detected QRS complex locations
    %   sensitivity: Detection sensitivity (if annotations are provided)
    %   positive_predictivity: Positive predictivity (if annotations are provided)

    % Extract the first channel of the ECG signal
    ecg_signal = signal(:, 1);

    % Step 1: Preprocessing - Apply bandpass filter (5-15 Hz)
    [b, a] = butter(5, [5 15]/(fs/2), 'bandpass');
    f = filtfilt(b, a, ecg_signal);

    % Define structuring element for morphological operations
    se_length = round(0.057 * fs);  % 57 ms (within the 55-60 ms range suggested)
    B = ones(se_length, 1);

    % Perform morphological operations
    f_erode = imerode(f, B);
    f_dilate = imdilate(f, B);
    f_open = imdilate(f_erode, B);
    f_close = imerode(f_dilate, B);

    % Detect peaks and valleys
    peaks = f - f_open;
    valleys = f - f_close;
    peaks_valleys = f - (f_open + f_close) / 2;

    % Find local maxima (potential QRS complexes)
    min_distance = round(0.2 * fs);  % Minimum 200 ms between peaks
    [~, qrs_peaks] = findpeaks(peaks_valleys, 'MinPeakDistance', min_distance);

    % Apply amplitude threshold to filter out noise
    threshold = 0.5 * max(peaks_valleys);
    qrs_peaks = qrs_peaks(peaks_valleys(qrs_peaks) > threshold);

    % Evaluate detection performance if annotations are provided
    if ~isempty(annotations)
        [sensitivity, positive_predictivity] = evaluate_detection(qrs_peaks, annotations, 0.15 * fs);
    else
        sensitivity = NaN;
        positive_predictivity = NaN;
        fprintf('Could not evaluate detection performance: annotations not available.\n');
    end
end