emg_agcl_all = readtable("Folder1-AgCl_EMG_properchair,stand,quadupndwn,squat,walkcrclandstrt_Default_none.xlsx.xlsx");
emg_hg_all = readtable("Folder1-PPHG EMG raw.xlsx");

emg_agcl = emg_agcl_all.AGCL_EMG_LEAD_I';
emg_agcl = emg_agcl(1:43850);

emg_hg = emg_hg_all.EMG_LEAD_I';
emg_hg = -emg_hg(1:43850);

[emg_hg, noise_segment_HG, SNR_dB_HG, f_agcl, P_clean_HG] = preprocess_emg_complete_adv(emg_hg', 200, true, true);


% Example variables (replace with your actual signals)
signal1 = emg_agcl; % first signal vector
signal2 = emg_hg; % second signal vector
fs = 200;      % sampling frequency in Hz

% Assuming signal1, signal2, and sampling frequency fs are defined

t1 = (0:length(signal1)-1) / fs;
t2 = (0:length(signal2)-1) / fs;

%% 1. Compute smooth envelopes (outline) using the Hilbert transform + smoothing

env1_raw = abs(hilbert(signal1));
env2_raw = abs(hilbert(signal2));

% Smooth the envelope using a moving average filter (adjust window size as needed)
smooth_window_sec = 0.03; % 50 ms smoothing window
smooth_window_samples = round(smooth_window_sec * fs);

env1 = movmean(env1_raw, smooth_window_samples);
env2 = movmean(env2_raw, smooth_window_samples);

%% 2. Plot signals with their smooth envelopes superimposed

figure('Name','Signals with Smooth Envelopes');

subplot(2,1,1);
plot(t1, signal1, 'b');
hold on;
plot(t1, env1, 'k', 'LineWidth', 1);
title('Signal 1 with Envelope Outline');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Signal 1', 'Envelope');
grid on;

subplot(2,1,2);
plot(t2, signal2, 'b');
hold on;
plot(t2, env2, 'k', 'LineWidth', 1);
title('Signal 2 with Envelope Outline');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Signal 2', 'Envelope');
grid on;

%% 3. Plot spectrograms side by side

window_length = round(0.5 * fs); % 0.5-second window
noverlap = round(0.4 * window_length);

figure('Name','Spectrogram Comparison');

subplot(1,2,1);
spectrogram(signal1, window_length, noverlap, 256, fs, 'yaxis');
title('Spectrogram of Signal 1');
ylim([0 fs/2]);
colorbar;

subplot(1,2,2);
spectrogram(signal2, window_length, noverlap, 256, fs, 'yaxis');
title('Spectrogram of Signal 2');
ylim([0 fs/2]);
colorbar;

%% 4. Plot scalograms (Continuous Wavelet Transform) side by side
figure('Name','Scalogram Comparison');

% Create subplot axes for first scalogram
ax1 = subplot(1,2,1);
cwt(signal1, fs, 'Parent', ax1);
title(ax1, 'Scalogram of Signal 1');

% Create subplot axes for second scalogram
ax2 = subplot(1,2,2);
cwt(signal2, fs, 'Parent', ax2);
title(ax2, 'Scalogram of Signal 2');


% --- Example usage ---

% Compute peak-to-peak values for both signals
ptp_signal1 = compute_peak_to_peak(signal1(24449:33615), fs);
ptp_signal2 = compute_peak_to_peak(signal2(26540:36071), fs);

% Calculate mean peak-to-peak voltage
mean_ptp1 = mean(ptp_signal1);
mean_ptp2 = mean(ptp_signal2);

% Bar plot comparison
figure;
bar([mean_ptp1, mean_ptp2]);
set(gca, 'XTickLabel', {'Signal 1', 'Signal 2'});
ylabel('Mean Peak-to-Peak Voltage');
title('Comparison of Mean Peak-to-Peak Voltage');
grid on;
