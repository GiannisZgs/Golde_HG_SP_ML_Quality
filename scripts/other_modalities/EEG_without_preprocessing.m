%% EEG Analysis Without Preprocessing Pipeline
% Analyzes raw EEG data without filtering or artifact removal
% Performs frequency domain analysis including PSD, spectrogram, and CWT
% Evaluates delta band power differences between conditions (eyes open vs. closed)

%% Combined EEG Analysis Script with PSD Comparison Plot

%% User Inputs
raw_eeg = amplitude_uV_HG(1:24000)'; % Replace with your actual continuous EEG data vector
fs = 200;                % Sampling frequency in Hz

% Define your condition time ranges in seconds
conditions = struct();
conditions.baseline = [0 10; 31 40]; % Baseline intervals
conditions.thinking = [11 30];                 % Thinking & eyes open
conditions.eyes_closed = [41 120];             % Eyes closed & rest

%% === Step 1: Exploratory Analysis ===

% Time vector for plotting
t = (0:length(raw_eeg)-1) / fs;


eeg_filtered_temp = raw_eeg;

% Plot raw and filtered EEG
figure('Name','Raw and Filtered EEG');
subplot(2,1,1);
plot(t, raw_eeg);
title('Raw EEG Signal');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(2,1,2);
plot(t, eeg_filtered_temp);
title('Filtered EEG Signal (1-40 Hz)');
xlabel('Time (s)'); ylabel('Amplitude');

% Compute and plot PSD for whole recording
window = hamming(floor(fs*2));
noverlap = floor(length(window)/2);
[pxx,freqs] = pwelch(eeg_filtered_temp, window, noverlap, 0:0.1:fs/2, fs);

figure('Name','PSD of Whole EEG Signal');
plot(freqs, 10*log10(pxx));
xlim([0 40]);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density - Whole EEG Signal');
grid on;

% Spectrogram visualization
figure('Name','EEG Spectrogram');
spectrogram(eeg_filtered_temp, window, noverlap, 256, fs, 'yaxis');
title('EEG Spectrogram');
ylim([1 40]);
colorbar;

% Continuous Wavelet Transform (scalogram)
figure('Name','EEG Scalogram (CWT)');
cwt(raw_eeg, fs);
title('EEG Scalogram (Continuous Wavelet Transform)-HG');
ylim([0 40]);

%% === Step 2: Alpha Power Epoching, Statistics, and Time-Resolved Alpha Power ===

% Use filtered EEG for analysis
eeg_filtered = eeg_filtered_temp;

% Epoching parameters
epoch_sec = 2;                  % epoch length in seconds
epoch_samples = epoch_sec * fs; % samples per epoch
overlap_sec = 1;                % overlap in seconds
step_samples = (epoch_sec - overlap_sec) * fs;

% Helper: convert seconds to samples
sec2samples = @(t) round(t*fs);

% Initialize alpha power arrays and PSD variables
alpha_baseline = [];
pxx_baseline = [];
f_baseline = [];

for i = 1:size(conditions.baseline,1)
    idx_start = sec2samples(conditions.baseline(i,1)) + 1;
    idx_end = sec2samples(conditions.baseline(i,2));
    segment = eeg_filtered(idx_start:idx_end);
    n_epochs = floor((length(segment) - epoch_samples) / step_samples) + 1;
    for ep = 1:n_epochs
        start_idx = (ep-1)*step_samples + 1;
        epoch_data = segment(start_idx : start_idx + epoch_samples - 1);
        [pxx_epoch,f_epoch] = pwelch(epoch_data, epoch_samples, 0, 0:0.1:fs/2, fs);
        alpha_band = (f_epoch >= 0.5) & (f_epoch <= 4);
        alpha_baseline(end+1,1) = trapz(f_epoch(alpha_band), pxx_epoch(alpha_band));
    end
    % For PSD comparison plot, compute PSD on entire baseline segment (once)
    if isempty(pxx_baseline)
        [pxx_baseline,f_baseline] = pwelch(segment, window, noverlap, 0:0.1:fs/2, fs);
    end
end

% Thinking (Eyes Open)
idx_start = sec2samples(conditions.thinking(1,1)) + 1;
idx_end = sec2samples(conditions.thinking(1,2));
segment = eeg_filtered(idx_start:idx_end);
alpha_thinking = [];
pxx_open = [];
f_open = [];
n_epochs = floor((length(segment) - epoch_samples) / step_samples) + 1;
for ep = 1:n_epochs
    start_idx = (ep-1)*step_samples + 1;
    epoch_data = segment(start_idx : start_idx + epoch_samples - 1);
    [pxx_epoch,f_epoch] = pwelch(epoch_data, epoch_samples, 0, 0:0.1:fs/2, fs);
    alpha_band = (f_epoch >= 0.5) & (f_epoch <= 4);
    alpha_thinking(end+1,1) = trapz(f_epoch(alpha_band), pxx_epoch(alpha_band));
end
[pxx_open,f_open] = pwelch(segment, window, noverlap, 0:0.1:fs/2, fs);

% Eyes Closed & Rest
idx_start = sec2samples(conditions.eyes_closed(1,1)) + 1;
idx_end = sec2samples(conditions.eyes_closed(1,2));
segment = eeg_filtered(idx_start:idx_end);
alpha_eyes_closed = [];
pxx_closed = [];
f_closed = [];
n_epochs = floor((length(segment) - epoch_samples) / step_samples) + 1;
for ep = 1:n_epochs
    start_idx = (ep-1)*step_samples + 1;
    epoch_data = segment(start_idx : start_idx + epoch_samples - 1);
    [pxx_epoch,f_epoch] = pwelch(epoch_data, epoch_samples, 0, 0:0.1:fs/2, fs);
    alpha_band = (f_epoch >= 0.5) & (f_epoch <= 4);
    alpha_eyes_closed(end+1,1) = trapz(f_epoch(alpha_band), pxx_epoch(alpha_band));
end
[pxx_closed,f_closed] = pwelch(segment, window, noverlap, 0:0.1:fs/2, fs);

%% Step 8: Plot PSD Comparison Between Epochs
figure('Name','PSD Comparison - Epochs');
plot(f_open, 10*log10(pxx_open), 'b', 'DisplayName', 'Thinking & Eyes Open'); hold on;
plot(f_closed, 10*log10(pxx_closed), 'r', 'DisplayName', 'Eyes Closed & Rest');
%plot(f_baseline, 10*log10(pxx_baseline), 'k', 'DisplayName', 'Baseline (No Experiment)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD Comparison Between Conditions');
xlim([0 40]);
legend('show');
grid on;

%% Statistical test: paired t-test between eyes closed and thinking
min_len = min(length(alpha_thinking), length(alpha_eyes_closed));
[h, p, ~, stats] = ttest(alpha_eyes_closed(1:min_len), alpha_thinking(1:min_len));
fprintf('Paired t-test result: p = %.4f, t = %.3f\n', p, stats.tstat);
if h == 1
    disp('Significant alpha power difference between Eyes Closed and Thinking (Eyes Open).');
else
    disp('No significant alpha power difference detected.');
end

% Boxplot of alpha power distributions
figure('Name','Delta Power Distribution by Condition');
boxplot([alpha_thinking(1:min_len), alpha_eyes_closed(1:min_len)], {'Thinking & Eyes Open', 'Eyes Closed & Rest'});
ylabel('Delta Power (0.5-4 Hz)');
title('Alpha Power Distribution by Condition');
grid on;

% Time-resolved alpha power with sliding windows
window_length_sec = 2;
step_size_sec = 0.5;
window_samples = round(window_length_sec * fs);
step_samples_full = round(step_size_sec * fs);

num_windows = floor((length(eeg_filtered) - window_samples) / step_samples_full) + 1;
alpha_power_time = zeros(1, num_windows);
time_vec = zeros(1, num_windows);

for i = 1:num_windows
    idx_start = (i-1)*step_samples_full + 1;
    idx_end = idx_start + window_samples - 1;
    segment = eeg_filtered(idx_start:idx_end);
    [pxx,f] = pwelch(segment, window_samples, 0, 0:0.1:fs/2, fs);
    alpha_band = (f >= 0.5) & (f <= 4);
    alpha_power_time(i) = trapz(f(alpha_band), pxx(alpha_band));
    time_vec(i) = (idx_start + idx_end)/(2*fs);
end

% Plot time-resolved alpha power with condition overlays
figure('Name','Time-Resolved Delta Power');
plot(time_vec, alpha_power_time, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Delta Power (0.5-4 Hz)');
title('Time-Resolved Delta Power-HG');
grid on;
hold on;

% Shade baseline intervals
y_limits = ylim;
for i = 1:size(conditions.baseline,1)
    x = conditions.baseline(i,:);
    fill([x(1) x(2) x(2) x(1)], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
% Shade thinking interval
x = conditions.thinking(1,:);
fill([x(1) x(2) x(2) x(1)], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
    [0.7 0.9 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% Shade eyes closed interval
x = conditions.eyes_closed(1,:);
fill([x(1) x(2) x(2) x(1)], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
    [0.7 0.7 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

legend('Delta Power', 'Baseline', 'Eyes Closed & Rest','Thinking & Eyes Open');
hold off;
