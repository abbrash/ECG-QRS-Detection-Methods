function [qrs_peaks, sensitivity, positive_predictivity] = detect_qrs_crossing(signal, Fs, ann)
    % QRS detection using zero crossing counts
    % Input:
    %   signal: input ECG signal
    %   Fs: sampling frequency in Hz
    %   ann: annotated QRS locations (optional)
    % Output:
    %   qrs_peaks: detected QRS complex locations
    %   sensitivity: detection sensitivity (if annotations provided)
    %   positive_predictivity: positive predictivity (if annotations provided)

    % Parameters
    window_size = round(0.15 * Fs);  % 150 ms window
    lambda_K = 0.98;  % forgetting factor for amplitude estimation
    lambda_D = 0.98;  % forgetting factor for feature signal
    lambda_theta = 0.98;  % forgetting factor for threshold
    c = 4;  % constant gain for high-frequency sequence

    % Band-pass filter (18-35 Hz)
    [b, a] = butter(3, [18 35] / (Fs/2), 'bandpass');
    ecg_filtered = filtfilt(b, a, signal);

    % Nonlinear transform
    y = abs(ecg_filtered);

    % Amplitude estimation and high-frequency sequence addition
    K = zeros(size(y));
    z = zeros(size(y));
    for n = 2:length(y)
        K(n) = lambda_K * K(n-1) + (1 - lambda_K) * abs(y(n));
        z(n) = y(n) + c * K(n) * (-1)^n;
    end

    % Zero crossing detection and counting
    zc = zeros(size(z));
    zc(2:end) = (z(1:end-1) .* z(2:end) <= 0);

    % Feature signal computation
    D = zeros(size(zc));
    for n = window_size:length(zc)
        D(n) = lambda_D * D(n-1) + (1 - lambda_D) * sum(zc(n-window_size+1:n));
    end

    % Adaptive thresholding
    theta = zeros(size(D));
    for n = 2:length(D)
        theta(n) = lambda_theta * theta(n-1) + (1 - lambda_theta) * D(n);
    end

    % Event detection
    events = D < theta;
    
    % Combine close events
    min_distance = round(0.2 * Fs);  % 200 ms minimum distance between QRS complexes
    events = combine_close_events(events, min_distance);

    % Finding QRS locations
    [~, qrs_peaks] = findpeaks(-D, 'MinPeakHeight', -max(theta), 'MinPeakDistance', min_distance);

    % Adjusting QRS locations to original signal
    for i = 1:length(qrs_peaks)
        [~, max_loc] = max(abs(signal(max(1, qrs_peaks(i)-round(0.1*Fs)):min(length(signal), qrs_peaks(i)+round(0.1*Fs)))));
        qrs_peaks(i) = qrs_peaks(i) + max_loc - round(0.1*Fs) - 1;
    end

    % Calculating performance metrics if annotations are provided
    if nargin > 2 && ~isempty(ann)
        [sensitivity, positive_predictivity] = calculate_performance(qrs_peaks, ann, Fs);
    else
        sensitivity = NaN;
        positive_predictivity = NaN;
    end
end

function combined_events = combine_close_events(events, min_distance)
    % Combine events that are closer than min_distance
    combined_events = events;
    event_starts = find(diff([0; events]) == 1);
    event_ends = find(diff([events; 0]) == -1);
    
    for i = 2:length(event_starts)
        if event_starts(i) - event_ends(i-1) < min_distance
            combined_events(event_ends(i-1):event_starts(i)) = 1;
        end
    end
end

function [sensitivity, positive_predictivity] = calculate_performance(detected_qrs, annotated_qrs, Fs)
    % Calculate performance metrics
    tolerance = round(0.15 * Fs);  % 150 ms tolerance window
    
    true_positives = 0;
    for i = 1:length(annotated_qrs)
        if any(abs(detected_qrs - annotated_qrs(i)) <= tolerance)
            true_positives = true_positives + 1;
        end
    end
    
    false_negatives = length(annotated_qrs) - true_positives;
    false_positives = length(detected_qrs) - true_positives;
    
    sensitivity = true_positives / (true_positives + false_negatives);
    positive_predictivity = true_positives / (true_positives + false_positives);
end