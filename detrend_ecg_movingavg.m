function ecg_detrended = detrend_ecg_movingavg(ecg_data, window_size)
    ecg_detrended = zeros(size(ecg_data));
    for ch = 1:size(ecg_data, 2)
        baseline = movmean(ecg_data(:, ch), window_size);
        ecg_detrended(:, ch) = ecg_data(:, ch) - baseline;
    end
end
