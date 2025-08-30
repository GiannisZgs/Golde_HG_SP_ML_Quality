function ecg_detrended = detrend_ecg_highpass(ecg_data, fs, cutoff, order)
    % Removes baseline wander from ECG signals using a high-pass Butterworth filter
    % Applies filtfilt for zero-phase filtering across all channels

    % High-pass filter parameters
    [b, a] = butter(order, cutoff / (fs / 2), 'high');
    
    % Apply filter to each channel
    ecg_detrended = zeros(size(ecg_data));
    for ch = 1:size(ecg_data, 2)
        ecg_detrended(:, ch) = filtfilt(b, a, ecg_data(:, ch));
    end
end
