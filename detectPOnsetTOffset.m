function [p_onset_idx, t_offset_idx] = detectPOnsetTOffset(ecg_signal, sampling_frequency, r_peaks,p_search_window_ms,t_search_window_ms)
%detectPOnsetTOffset Detects P-wave onset and T-wave offset in an ECG signal.
%   This function assumes R-peaks have already been detected and uses them
%   as reference points.
%
%   Inputs:
%     ecg_signal: A 1D array containing the ECG signal.
%     sampling_frequency: The sampling frequency of the ECG signal in Hz.
%     r_peaks: A 1D array of indices corresponding to the detected R-peak locations.
%
%   Outputs:
%     p_onset_idx: Indices of detected P-wave onset points.
%     t_offset_idx: Indices of detected T-wave offset points.

% --- Parameters (Adjust these for optimal performance based on your data) ---
% Bandpass filter for P and T waves
f_low_pt = 0.5;   % Hz, lower cutoff frequency for P/T waves
f_high_pt = 20;   % Hz, upper cutoff frequency for P/T waves (P/T waves are typically lower freq)

% Thresholding parameters for slope-based detection
% These are crucial for identifying the start/end of the waves based on slope changes.
% Values are typically a fraction of the maximum derivative within the wave.
p_onset_slope_factor = 0.05; % A smaller factor means earlier onset detection (more sensitive)
t_offset_slope_factor = 0.05; % A smaller factor means earlier offset detection (more sensitive)

% --- Preprocessing: Bandpass Filter for P and T waves ---
% Apply a Butterworth bandpass filter to isolate the frequency components
% relevant to P and T waves, while reducing baseline wander and high-freq noise.
Wn_pt = [f_low_pt, f_high_pt] / (sampling_frequency / 2);
[b_pt, a_pt] = butter(4, Wn_pt, 'bandpass'); % 4th order filter, good balance of steepness and ripple
filtered_ecg_pt = filtfilt(b_pt, a_pt, ecg_signal); % Zero-phase filtering

% Calculate the derivative of the filtered signal.
% The derivative helps highlight slope changes, which are key for onset/offset.
deriv_filtered_ecg_pt = diff([0; filtered_ecg_pt(:)]); % Add a zero for length matching

% --- Initialize Output Arrays ---
p_onset_idx = NaN(1, length(r_peaks)); % Pre-allocate with NaN in case a wave isn't found
t_offset_idx = NaN(1, length(r_peaks));

% --- Loop Through Each Detected R-peak to Find P and T Waves ---
for i = 1:length(r_peaks)
    r_idx = r_peaks(i); % Current R-peak index

    % --- P-wave Onset Detection ---
    % Define the search window before the R-peak.
    % The P-wave search stops a little before the R-peak to avoid QRS interference.
    p_segment_end = r_idx - round(0.02 * sampling_frequency); % End search 20ms before R-peak
    p_segment_start = max(1, r_idx - round(p_search_window_ms / 1000 * sampling_frequency));

    % Ensure a valid search window exists
    if p_segment_start < p_segment_end
        current_p_segment = filtered_ecg_pt(p_segment_start:p_segment_end);
        current_p_segment_deriv = deriv_filtered_ecg_pt(p_segment_start:p_segment_end);

        % Find the P-wave peak within the search window.
        % We look for the largest local maximum (P-peak).
        [p_peaks_val, p_peaks_loc_rel] = findpeaks(current_p_segment, 'MinPeakHeight', 0.01 * max(abs(current_p_segment)), 'NPeaks', 1, 'SortStr','descend');
        
        if ~isempty(p_peaks_val)
            p_peak_idx_global = p_segment_start + p_peaks_loc_rel(1) - 1; % Global index of P-peak
            
            % Search left from P-peak for onset (where the slope begins to increase)
            % The onset is typically when the derivative crosses a small threshold
            % relative to the maximum derivative of the P wave.
            search_from_p_peak_left = p_peak_idx_global:-1:p_segment_start;
            max_deriv_p_wave = max(deriv_filtered_ecg_pt(search_from_p_peak_left));
            
            found_p_onset = false;
            for k = search_from_p_peak_left
                if deriv_filtered_ecg_pt(k) > (p_onset_slope_factor * max_deriv_p_wave)
                    % The point where the positive slope begins to be significant.
                    p_onset_idx(i) = k;
                    found_p_onset = true;
                    break;
                end
            end
            
            if ~found_p_onset && ~isempty(search_from_p_peak_left)
                p_onset_idx(i) = search_from_p_peak_left(end); % Fallback to start of search if no clear onset
            end
        end
    end

    % --- T-wave Offset Detection ---
    % Define the search window after the R-peak (or QRS offset if available).
    % Start slightly after the QRS complex to avoid interference.
    t_segment_start = min(length(ecg_signal), r_idx + round(0.05 * sampling_frequency)); % Start 50ms after R
    t_segment_end = min(length(ecg_signal), r_idx + round(t_search_window_ms / 1000 * sampling_frequency));

    % Ensure a valid search window exists
    if t_segment_start < t_segment_end
        current_t_segment = filtered_ecg_pt(t_segment_start:t_segment_end);
        current_t_segment_deriv = deriv_filtered_ecg_pt(t_segment_start:t_segment_end);

        % Find the T-wave peak within the search window.
        % We look for the largest local maximum (T-peak) if T-wave is positive.
        [t_peaks_val, t_peaks_loc_rel] = findpeaks(current_t_segment, 'MinPeakHeight', 0.1 * max(abs(current_t_segment)), 'NPeaks', 1, 'SortStr','descend');
        
        if ~isempty(t_peaks_val)
            t_peak_idx_global = t_segment_start + t_peaks_loc_rel(1) - 1; % Global index of T-peak
            
            % Search right from T-peak for offset (where the slope flattens to baseline)
            % The offset is typically when the derivative approaches zero or crosses a small threshold.
            search_from_t_peak_right = t_peak_idx_global:1:t_segment_end;
            min_deriv_t_wave = min(deriv_filtered_ecg_pt(search_from_t_peak_right)); % Max negative slope of T wave
            
            found_t_offset = false;
            for k = search_from_t_peak_right
                % We're looking for where the negative slope of the T-wave returns to flat
                if abs(deriv_filtered_ecg_pt(k)) < abs(t_offset_slope_factor * min_deriv_t_wave)
                    t_offset_idx(i) = k;
                    found_t_offset = true;
                    break;
                end
            end
            
            if ~found_t_offset && ~isempty(search_from_t_peak_right)
                t_offset_idx(i) = search_from_t_peak_right(end); % Fallback to end of search if no clear offset
            end
        end
    end
end

% Remove any NaN values (beats where P/T wave couldn't be definitively found)
p_onset_idx = p_onset_idx(~isnan(p_onset_idx));
t_offset_idx = t_offset_idx(~isnan(t_offset_idx));

end