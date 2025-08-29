function ecg_detrended = detrend_ecg_highpass(ecg_data, fs, cutoff, order)

    % High-pass filter parameters
    [b, a] = butter(order, cutoff / (fs / 2), 'high');
    
    % Apply filter to each channel
    ecg_detrended = zeros(size(ecg_data));
    for ch = 1:size(ecg_data, 2)
        ecg_detrended(:, ch) = filtfilt(b, a, ecg_data(:, ch));
    end
end
