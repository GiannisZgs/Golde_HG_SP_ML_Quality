%% Full EOG Preprocessing + Peak Analysis + Statistics + Plots

% User inputs
raw_eog = amplitude_uV(7111:10289)'; % Replace with your actual raw EOG vector
fs = 200; % Sampling frequency in Hz

t = (0:length(raw_eog)-1)/fs; % time vector

%% 1. Bandpass filter (0.5 - 30 Hz)
[b_bp, a_bp] = butter(4, [0.5 30]/(fs/2), 'bandpass');
eog_bp = filtfilt(b_bp, a_bp, raw_eog);

%% 2. Notch filter at 50 Hz
f0 = 50;
Q = 30;
bw = f0 / Q;
[b_notch, a_notch] = iirnotch(f0/(fs/2), bw/(fs/2));
eog_notch = filtfilt(b_notch, a_notch, eog_bp);

%% 3. Baseline correction (mean removal)
eog_basecorr = eog_notch - mean(eog_notch);

%% 4. Detrend signal
eog_detrended = detrend(eog_basecorr);

%% 5. Wavelet denoising
eog_denoised = wdenoise(eog_detrended, 3, 'Wavelet', 'db4', ...
    'DenoisingMethod', 'Bayes', 'ThresholdRule', 'Soft', 'NoiseEstimate', 'LevelDependent');

%% 6. Median filtering
median_filter_window = round(0.03 * fs); % 50 ms window
eog_median = medfilt1(eog_denoised, median_filter_window);

% Plot preprocessing stages
figure;
subplot(4,1,1); plot(t, raw_eog); title('Raw EOG Signal'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,1,2); plot(t, eog_basecorr); title('Filtered + Notch + Baseline Corrected EOG'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,1,3); plot(t, eog_denoised); title('After Wavelet Denoising'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(4,1,4); plot(t, eog_median); title('After Median Filtering'); xlabel('Time (s)'); ylabel('Amplitude');

eog_clean = eog_median; % final cleaned EOG

%% Define activity intervals (seconds)
activity_times.look_up = [0, 5];
activity_times.blink_3 = [5.2, 7.3];
activity_times.look_down = [9.7, 15.5];

% Epoch parameters
epoch_sec = 1; % seconds
epoch_samples = epoch_sec * fs;
overlap_sec = 0.5;
step_samples = (epoch_sec - overlap_sec) * fs;

sec2samples = @(t) round(t*fs);

%% Function to compute peaks per epoch
% function [max_pks, min_pks, ptp_vals] = compute_peaks(signal_segment, epoch_samples, step_samples)
%     n_epochs = floor((length(signal_segment) - epoch_samples)/step_samples) + 1;
%     max_pks = zeros(n_epochs,1);
%     min_pks = zeros(n_epochs,1);
%     ptp_vals = zeros(n_epochs,1);
%     for i = 1:n_epochs
%         idx_start = (i-1)*step_samples + 1;
%         idx_end = idx_start + epoch_samples - 1;
%         epoch = signal_segment(idx_start:idx_end);
%         max_pks(i) = max(epoch);
%         min_pks(i) = min(epoch);
%         ptp_vals(i) = max_pks(i) - min_pks(i);
%     end
% end

% Extract look up and look down segments
seg_up = eog_clean(sec2samples(activity_times.look_up(1))+1 : sec2samples(activity_times.look_up(2)));
seg_down = eog_clean(sec2samples(activity_times.look_down(1))+1 : sec2samples(activity_times.look_down(2)));

[max_up, min_up, ptp_up] = compute_peaks(seg_up, epoch_samples, step_samples);
[max_down, min_down, ptp_down] = compute_peaks(seg_down, epoch_samples, step_samples);

%% Blink detection and PtP in blink segment

seg_blink = eog_clean(sec2samples(activity_times.blink_3(1))+1 : sec2samples(activity_times.blink_3(2)));

min_peak_height = 0.3 * max(abs(seg_blink));
min_peak_distance = round(0.2 * fs);

[pos_peaks, pos_locs] = findpeaks(seg_blink, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);
[neg_peaks, neg_locs] = findpeaks(-seg_blink, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);

pos_locs = pos_locs(:);
neg_locs = neg_locs(:);
all_locs = sort([pos_locs; neg_locs]);


win_len = round(0.4 * fs);
half_win = round(win_len / 2);

ptp_blinks = zeros(length(all_locs), 1);

for i = 1:length(all_locs)
    idx_center = all_locs(i);
    idx_start = max(idx_center - half_win, 1);
    idx_end = min(idx_center + half_win, length(seg_blink));
    blink_waveform = seg_blink(idx_start:idx_end);
    ptp_blinks(i) = max(blink_waveform) - min(blink_waveform);
end

% Use mean peak values for plotting
max_blink = mean(pos_peaks);
min_blink = mean(-neg_peaks);
ptp_blink = mean(ptp_blinks);

%% Prepare data for bar plots

max_all = [mean(max_up), max_blink, mean(max_down)];
min_all = [mean(min_up), min_blink, mean(min_down)];
max_all = max_all(:)';
min_all = min_all(:)';

ptp_all = [mean(ptp_up), ptp_blink, mean(ptp_down)];
ptp_all = ptp_all(:)';

labels = {'Look Up', 'Blink 3x', 'Look Down'};

% Bar plot: mean positive and absolute mean negative peaks
figure;
bar([max_all; abs(min_all)]');
set(gca, 'XTickLabel', labels);
legend('Mean Max Peak', 'Mean |Min Peak|');
ylabel('Amplitude');
title('Mean Positive and Negative Peak Amplitudes per Activity');
grid on;

% Bar plot: mean peak-to-peak voltages
figure;
bar(ptp_all);
set(gca, 'XTickLabel', labels);
ylabel('Mean Peak-to-Peak Voltage');
title('Mean Peak-to-Peak Voltage per Activity');
grid on;

%% Statistical Analysis: One-way ANOVA + Tukey post-hoc

% Combine all data for ANOVA
all_ptp = [ptp_up(:); ptp_blinks(:); ptp_down(:)];
group_ptp = [ ...
    repmat({'Look Up'}, length(ptp_up), 1); ...
    repmat({'Blink 3x'}, length(ptp_blinks), 1); ...
    repmat({'Look Down'}, length(ptp_down), 1)];

all_max = [max_up(:); pos_peaks(:); max_down(:)];
group_max = [ ...
    repmat({'Look Up'}, length(max_up), 1); ...
    repmat({'Blink 3x'}, length(pos_peaks), 1); ...
    repmat({'Look Down'}, length(max_down), 1)];

all_min = [min_up(:); -neg_peaks(:); min_down(:)];
group_min = [ ...
    repmat({'Look Up'}, length(min_up), 1); ...
    repmat({'Blink 3x'}, length(neg_peaks), 1); ...
    repmat({'Look Down'}, length(min_down), 1)];


% ANOVA for peak-to-peak
[p_ptp, ~, stats_ptp] = anova1(all_ptp, group_ptp);
fprintf('ANOVA peak-to-peak p-value: %.4f\n', p_ptp);
if p_ptp < 0.05
    figure; multcompare(stats_ptp, 'CType', 'tukey-kramer');
end

% ANOVA for max peak
[p_max, ~, stats_max] = anova1(all_max, group_max);
fprintf('ANOVA max peak p-value: %.4f\n', p_max);
if p_max < 0.05
    figure; multcompare(stats_max, 'CType', 'tukey-kramer');
end

% ANOVA for min peak
[p_min, ~, stats_min] = anova1(all_min, group_min);
fprintf('ANOVA min peak p-value: %.4f\n', p_min);
if p_min < 0.05
    figure; multcompare(stats_min, 'CType', 'tukey-kramer');
end

close all;