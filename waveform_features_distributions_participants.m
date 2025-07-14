clear
close all

manually_cleaned = 1;
if manually_cleaned
    data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\manually_cleaned_participant_id_results");
else
    data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\participant_id_results");
end
features = data.features;

sensors = ["AgCl","HG"];
fs = 200;

%% ECG waveform features Boxplots distributions
participants = fields(features);
peak_to_peak = cell(6,1);
qrs_duration = cell(6,1);
t_wave_amp = cell(6,1);
p_wave_amp = cell(6,1);

group_names = {'AgCl Channel 1', 'HG Channel 1', 'AgCl Channel 2', 'HG Channel 2', 'AgCl Channel 3', 'HG Channel 3'};  
for p = 1:length(participants)
    participant = participants(p); participant = participant{1};

    for s = 1:2
        if isfield(features.(participant),(sensors(s)))
            for ch = 1:3
                if isfield(features.(participant).(sensors(s)),(['ch',num2str(ch)]))
                    %% Peak-to-peak amplitude
                    p2p_val = features.(participant).(sensors(s)).(['ch',num2str(ch)]).peak_to_peak;
                    if ~isnan(p2p_val)
                        if s == 1 && ch == 1
                            peak_to_peak(1,1) = {[peak_to_peak{1,1} p2p_val]}; 
                        elseif s == 2 && ch == 1
                            peak_to_peak(2,1) = {[peak_to_peak{2,1} p2p_val]}; 
                        elseif s == 1 && ch == 2
                            peak_to_peak(3,1) = {[peak_to_peak{3,1} p2p_val]}; 
                        elseif s == 2 && ch == 2
                            peak_to_peak(4,1) = {[peak_to_peak{4,1} p2p_val]}; 
                        elseif s == 1 && ch == 3
                            peak_to_peak(5,1) = {[peak_to_peak{5,1} p2p_val]}; 
                        elseif s == 2 && ch == 3
                            peak_to_peak(6,1) = {[peak_to_peak{6,1} p2p_val]}; 
                        end
                    end
                    %% QRS Duration
                    qrs_dur_val = 1000*features.(participant).(sensors(s)).(['ch',num2str(ch)]).qrs_duration;
                    if ~isnan(qrs_dur_val)
                        if s == 1 && ch == 1
                            qrs_duration(1,1) = {[qrs_duration{1,1} qrs_dur_val]}; 
                        elseif s == 2 && ch == 1
                            qrs_duration(2,1) = {[qrs_duration{2,1} qrs_dur_val]}; 
                        elseif s == 1 && ch == 2
                            qrs_duration(3,1) = {[qrs_duration{3,1} qrs_dur_val]}; 
                        elseif s == 2 && ch == 2
                            qrs_duration(4,1) = {[qrs_duration{4,1} qrs_dur_val]}; 
                        elseif s == 1 && ch == 3
                            qrs_duration(5,1) = {[qrs_duration{5,1} qrs_dur_val]}; 
                        elseif s == 2 && ch == 3
                            qrs_duration(6,1) = {[qrs_duration{6,1} qrs_dur_val]}; 
                        end
                    end
                    %% T-wave Amplitude
                    t_val = features.(participant).(sensors(s)).(['ch',num2str(ch)]).t_wave_amp;
                    if ~isnan(t_val)
                        if s == 1 && ch == 1
                            t_wave_amp(1,1) = {[t_wave_amp{1,1} t_val]}; 
                        elseif s == 2 && ch == 1
                            t_wave_amp(2,1) = {[t_wave_amp{2,1} t_val]}; 
                        elseif s == 1 && ch == 2
                            t_wave_amp(3,1) = {[t_wave_amp{3,1} t_val]}; 
                        elseif s == 2 && ch == 2
                            t_wave_amp(4,1) = {[t_wave_amp{4,1} t_val]}; 
                        elseif s == 1 && ch == 3
                            t_wave_amp(5,1) = {[t_wave_amp{5,1} t_val]}; 
                        elseif s == 2 && ch == 3
                            t_wave_amp(6,1) = {[t_wave_amp{6,1} t_val]}; 
                        end
                    end
                    %% P-wave Amplitude
                    p_val = features.(participant).(sensors(s)).(['ch',num2str(ch)]).p_wave_amp;
                    if ~isnan(p_val)
                        if s == 1 && ch == 1
                            p_wave_amp(1,1) = {[p_wave_amp{1,1} p_val]}; 
                        elseif s == 2 && ch == 1
                            p_wave_amp(2,1) = {[p_wave_amp{2,1} p_val]}; 
                        elseif s == 1 && ch == 2
                            p_wave_amp(3,1) = {[p_wave_amp{3,1} p_val]}; 
                        elseif s == 2 && ch == 2
                            p_wave_amp(4,1) = {[p_wave_amp{4,1} p_val]}; 
                        elseif s == 1 && ch == 3
                            p_wave_amp(5,1) = {[p_wave_amp{5,1} p_val]}; 
                        elseif s == 2 && ch == 3
                            p_wave_amp(6,1) = {[p_wave_amp{6,1} p_val]}; 
                        end
                    end
                end
            end
        end
    end
end

%% Peak-to-peak Group Boxplot
nGroups = length(peak_to_peak);
peak2peak_p_vals = NaN(nGroups);

% Perform pairwise Wilcoxon rank-sum tests (non-parametric)
for i = 1:nGroups
    for j = i+1:nGroups
        if j-i == 1 && mod(i,2) == 1
            peak2peak_p_vals(i,j) = ranksum(peak_to_peak{i}, peak_to_peak{j});
        end
    end
end

peak2peak_labels = cellfun(@(x, i) repmat(i, numel(x), 1), peak_to_peak, num2cell(1:numel(peak_to_peak))', 'UniformOutput', false);
peak2peak_labels = vertcat(peak2peak_labels{:})';
peak2peak = cell2mat(peak_to_peak');
figure;
boxplot(peak2peak, peak2peak_labels, 'Labels', group_names)
title('Peak-to-Peak ECG Profile Amplitude - Distribution of Participants');
ylabel('Peak-to-Peak Amplitude (mV)');
hold on;
group_positions = 1:nGroups;
y_offset = 0.05 * range(ylim); 
p_significance_on_plot(peak2peak_p_vals, y_offset, group_positions);

%% QRS Duration Group Boxplot
nGroups = length(qrs_duration);
qrs_dur_p_vals = NaN(nGroups);

% Perform pairwise Wilcoxon rank-sum tests (non-parametric)
for i = 1:nGroups
    for j = i+1:nGroups
        if j-i == 1 && mod(i,2) == 1
            qrs_dur_p_vals(i,j) = ranksum(qrs_duration{i}, qrs_duration{j});
        end
    end
end

qrs_labels = cellfun(@(x, i) repmat(i, numel(x), 1), qrs_duration, num2cell(1:numel(qrs_duration))', 'UniformOutput', false);
qrs_labels = vertcat(qrs_labels{:})';
qrs_duration = cell2mat(qrs_duration');
figure;
boxplot(qrs_duration, qrs_labels, 'Labels', group_names)
title('QRS Duration ECG Profile - Distribution of Participants');
ylabel('QRS Duration (ms)');
hold on;
group_positions = 1:nGroups;
y_offset = 0.05 * range(ylim); 
p_significance_on_plot(qrs_dur_p_vals, y_offset, group_positions);

%% T-wave Amplitude Group Boxplot
nGroups = length(t_wave_amp);
t_wave_p_vals = NaN(nGroups);

% Perform pairwise Wilcoxon rank-sum tests (non-parametric)
for i = 1:nGroups
    for j = i+1:nGroups
        if j-i == 1 && mod(i,2) == 1
            t_wave_p_vals(i,j) = ranksum(t_wave_amp{i}, t_wave_amp{j});
        end
    end
end

t_labels = cellfun(@(x, i) repmat(i, numel(x), 1), t_wave_amp, num2cell(1:numel(t_wave_amp))', 'UniformOutput', false);
t_labels = vertcat(t_labels{:})';
t_wave_amp = cell2mat(t_wave_amp');
figure;
boxplot(t_wave_amp, t_labels, 'Labels', group_names)
title('T-wave ECG Profile Amplitude - Distribution of Participants');
ylabel('Amplitude (mV)');
hold on;
group_positions = 1:nGroups;
y_offset = 0.05 * range(ylim); 
p_significance_on_plot(t_wave_p_vals, y_offset, group_positions);

%% P-wave Amplitude Group Boxplot
nGroups = length(p_wave_amp);
p_wave_p_vals = NaN(nGroups);

% Perform pairwise Wilcoxon rank-sum tests (non-parametric)
for i = 1:nGroups
    for j = i+1:nGroups
        if j-i == 1 && mod(i,2) == 1
            p_wave_p_vals(i,j) = ranksum(p_wave_amp{i}, p_wave_amp{j});
        end
    end
end

p_labels = cellfun(@(x, i) repmat(i, numel(x), 1), p_wave_amp, num2cell(1:numel(p_wave_amp))', 'UniformOutput', false);
p_labels = vertcat(p_labels{:})';
p_wave_amp = cell2mat(p_wave_amp');
figure;
boxplot(p_wave_amp, p_labels, 'Labels', group_names)
title('P-wave ECG Profile Amplitude - Distribution of Participants');
ylabel('Amplitude (mV)');
hold on;
group_positions = 1:nGroups;
y_offset = 0.05 * range(ylim); 
p_significance_on_plot(p_wave_p_vals, y_offset, group_positions);