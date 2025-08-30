%% ECG Data Extraction From Manually Cleaned Recordings
% Extracts and processes manually cleaned ECG data for participants 1, 5, 10, and 39
% Applies various signal processing filters and quality assessment metrics
% Compares AgCl and hydrogel electrodes for benchmarking performance

close all
clear;
%% Setup environment
[scripts_dir, ~, ~] = fileparts(pwd);
[root_dir, ~, ~] = fileparts(scripts_dir);
setup_script = fullfile(root_dir,'utils','setup_environment.m');
run(setup_script);

%% Path to the dataset goes here
data_path = 'C:\Users\giann\OneDrive\Desktop\ECG HG paper\Motion_artifact_cleaned_signals\manually_cleaned_ECG_HG_quality_dataset_no_filters.mat';

%% Parameters
save_files = true;
data_struct.fs_ecg = 200;
data_struct.fs_imp = 1;
window_size = data_struct.fs_ecg/data_struct.fs_imp;
num_windows = 2;
hp_cutoff = 0.05; %BW HP filter cutoff (Hz) for detrending
hp_order = 4; %BW HP filter order 
lp_cutoff = 90; %BW LP filter cutoff (Hz) - 100 is the Nyquist frequency
lp_order = 4;
notch_freq1 = 50; %Central bandstop frequency (UAE power line frequency)
notch_freq2 = 80; %Power line 1st harmonic artifact 
notch_bw = 2; 
thres_Imp = 100; %impedance > 100 kΩ considered not acceptable
use_filters = 1; %boolean, whether to use any filters in the analysis

%% Only for the below participants whose recordings were manually examined and cleaned
participants = ["1","5","10","39"];
signals = ["AgCl","HG1","HG2"];

data_struct = load(data_path); 
data_struct.fs_ecg = data_struct.data_struct.fs_ecg;
fieldss = fields(data_struct);

empty_fun = @(s) all(structfun(@isempty,s)); %checks if structure contains anything
no_AgCl = [];
no_HG1 = [];
no_HG2 = [];
no_HG = [];
no_all = [];
for p = 1:length(participants)
    p_struct = struct();
    flag_empty = zeros(length(signals),1);
    for s = 1:length(signals)
        if s==1
            signal_label = 'AgCl';
        elseif s==2
            signal_label = 'HG1';
        elseif s==3
            signal_label = 'HG2';
        end
        
        for f = 1:length(fieldss)
            fieldd = fieldss{f};
            if strcmp(fieldd,(['p' char(participants(p)) '_' signal_label '_ECG_ch3_clipped']))
                ECG_ch3_clipped = data_struct.(fieldd);
                ECG = data_struct.data_struct.(['p' char(participants(p))]).(signal_label).ECG;
            
                %% Align the other 2 channels using the timestamps of ch3
                %Time 0 is the first sample
                sample_dur = 1/data_struct.fs_ecg;
                start_time = 0; %Sample #1 is time 0
                try
                    sample_idxs = int64(time2num(ECG_ch3_clipped.Time) / sample_dur + 1);
                catch
                    sample_idxs = int64(time2num(ECG_ch3_clipped) / sample_dur + 1);
                end
                ECG = ECG(sample_idxs,:);
            else 
                continue
            end
            
            %% Rest of the pipeline remains the same except Impedance check
            if use_filters
                %ECG = ECG-mean(ECG);
                [blp, alp] = butter(lp_order, lp_cutoff / (data_struct.fs_ecg/2), 'low');
                ECG_lp = filtfilt(blp, alp, ECG);
                Wn1 = [notch_freq1 - notch_bw, notch_freq1 + notch_bw] / (data_struct.fs_ecg/2);
                Wn2 = [notch_freq2 - notch_bw, notch_freq2 + notch_bw] / (data_struct.fs_ecg/2);
                [bnotch1, anotch1] = butter(2, Wn1, 'stop');
                [bnotch2, anotch2] = butter(2, Wn2, 'stop');
                ECG_notch1 = filtfilt(bnotch1, anotch1, ECG_lp);
                ECG_just_notch1 = filtfilt(bnotch1, anotch1, ECG);
                ECG_notch2 = filtfilt(bnotch2, anotch2, ECG_notch1);
                ECG_just_notch2 = filtfilt(bnotch2, anotch2, ECG);
                ECG_just_hp_VLF = detrend_ecg_highpass(ECG,data_struct.fs_ecg,hp_cutoff,hp_order);
                ECG_filtered = detrend_ecg_movingavg(ECG_notch2,data_struct.fs_ecg); %1-second window size
                ECG_just_MA = detrend_ecg_movingavg(ECG,data_struct.fs_ecg); %
            else
                ECG_filtered = ECG;
            end

            %% Measure metrics to quantify the analogies of noise in AgCl and HG
            plot_flag = false;
            if use_filters 
                trend_deviation = deviation_from_noise_psd(ECG, ECG_just_MA, data_struct.fs_ecg, 'trend', plot_flag,signal_label,['p',char(participants(p))],true); 
                metrics.(['p',char(participants(p))]).(signal_label).trend_deviation = trend_deviation.trend_reduction_ratio;
                vlf_deviation = deviation_from_noise_psd(ECG, ECG_just_hp_VLF, data_struct.fs_ecg, 'vlf', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).vlf_deviation = vlf_deviation.vlf_reduction_ratio;
                low_freq_deviation = deviation_from_noise_psd(ECG, ECG_just_hp_VLF, data_struct.fs_ecg, 'low', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).low_freq_deviation = low_freq_deviation.low_reduction_ratio;
                mid_freq_deviation = deviation_from_noise_psd(ECG, ECG_just_hp_VLF, data_struct.fs_ecg, 'mid', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).mid_freq_deviation = mid_freq_deviation.mid_reduction_ratio;
                high_freq_deviation = deviation_from_noise_psd(ECG, ECG_just_hp_VLF, data_struct.fs_ecg, 'high', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).high_freq_deviation = high_freq_deviation.high_reduction_ratio;
                powerline_deviation = deviation_from_noise_psd(ECG, ECG_just_notch1, data_struct.fs_ecg, 'powerline', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).powerline_deviation = powerline_deviation.powerline_reduction_ratio;
                powerline_harmonic_deviation = deviation_from_noise_psd(ECG, ECG_just_notch2, data_struct.fs_ecg, 'powerline_harmonic', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).powerline_harmonic_deviation = powerline_harmonic_deviation.powerline_harmonic_reduction_ratio;
                ultra_high1_freq_deviation = deviation_from_noise_psd(ECG, ECG_lp, data_struct.fs_ecg, 'ultra_high1', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).ultra_high1_freq_deviation = ultra_high1_freq_deviation.ultra_high1_reduction_ratio;
                ultra_high2_freq_deviation = deviation_from_noise_psd(ECG, ECG_lp, data_struct.fs_ecg, 'ultra_high2', plot_flag,signal_label,['p',char(participants(p))],true);
                metrics.(['p',char(participants(p))]).(signal_label).ultra_high2_freq_deviation = ultra_high2_freq_deviation.ultra_high2_reduction_ratio;
            end
            plot_flag = true;
            if plot_flag && use_filters
                deviation_from_noise_psd(ECG, ECG_filtered, data_struct.fs_ecg, '', plot_flag,signal_label,['p',char(participants(p))],true); 
                close all;
            end

            if s==1
                raw_AgCl = ECG;
                proc_AgCl = ECG_filtered;
            elseif s==2
                raw_HG1 = ECG;
                proc_HG1 = ECG_filtered;
            elseif s==3
                raw_HG2 = ECG;
                proc_HG2 = ECG_filtered;
            end

            %% Apply a QRS detection algorithm on Channel 1
            try
                [QRS_indexes,filt_dat,int_dat,thF1,thI1] = pantompkins_qrs(ECG_filtered(:,1),data_struct.fs_ecg);
            catch ME
                if strcmp(ME.identifier,'MATLAB:unassignedOutputs')
                    QRS_indexes = [];
                    fprintf("\n This recording will be flagged for removal - QRS detection failed because of no signal in part or whole of the recording \n")
                else
                    rethrow(ME)
                end

            end
            %plot_pantompkins_qrs_detection(ECG_filtered(:,1),filt_dat,int_dat,QRS_indexes,thF1,thI1, data_struct.fs_ecg);

            %% Apply an impedance threshold criterion to discard ECG segments of bad quality
            ECG_good_qual = struct();
            imp_good_qual = struct();
            QRS_annot_good_qual = struct();
            for ch = 1:size(ECG,2)
                inds_start = (1:num_windows*window_size:size(ECG,1));
                inds_end = inds_start + relative_window_size*data_struct.fs_ecg - 1;

                %Find qrs indexes that lie within the [inds_start,inds_end] interval

                ecg_good_qual = [];
                qrs_annot_good_qual = cell(1,length(inds_start));
                breakpoints = [];
                inds_to_remove = 0;
                for i = 1:length(inds_start)
                    win_start = inds_start(i); win_end = inds_end(i);
                    if win_end - 1 > size(ECG,1)
                        inds_to_remove = inds_to_remove + 1;
                        continue
                    end
                    ecg_good_qual = [ecg_good_qual ECG_filtered(win_start:win_end,ch)];
                    breakpoints = [breakpoints [win_start;win_end]];
                    %Find qrs indexes that lie within the [inds_start,inds_end] interval    
                    qrs_higher = find(QRS_indexes>=win_start);
                    qrs_lower = find(QRS_indexes<=win_end);
                    qrs_common = intersect(qrs_higher,qrs_lower);
                    qrs_current = QRS_indexes(qrs_common);
                    if any(qrs_current > num_windows*window_size)
                        %Either all or none
                        qrs_current = qrs_current - win_start + 1;
                    end
                    qrs_annot_good_qual{i} = qrs_current; %+1
                end
                ECG_good_qual.(['ch',int2str(ch)]) = ecg_good_qual;
                QRS_annot_good_qual.(['ch',int2str(ch)]) = qrs_annot_good_qual(1:end-inds_to_remove);
                ECG_good_qual.(['ch',int2str(ch),'breakpoints']) = breakpoints;       
            end


            p_struct.(signal_label).ECG = ECG;
            p_struct.(signal_label).ECG_good_qual = ECG_good_qual;
            p_struct.(signal_label).QRS_annot_good_qual = QRS_annot_good_qual;
        end
    end
    plot_flag = true;
    if plot_flag && use_filters
        min_length = min([size(raw_AgCl,1), size(raw_HG1,1), size(raw_HG2,1)]);
        deviation_from_noise_psd(raw_AgCl(1:min_length,:), raw_HG1(1:min_length,:), data_struct.fs_ecg, '', plot_flag,'AgCl-HG1',['p',char(participants(p))],false); 
        deviation_from_noise_psd(raw_AgCl(1:min_length,:), raw_HG2(1:min_length,:), data_struct.fs_ecg, '', plot_flag,'AgCl-HG2',['p',char(participants(p))],false); 
        close all;
    end

    %% Concatenate the hydrogel signals 
    if ~flag_empty(2) || ~flag_empty(3)
        if ~flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch1))) && ~sum(any(isnan(p_struct.HG2.ECG_good_qual.ch1))))
            ch1_concat_ecg = p_struct.HG2.ECG_good_qual.ch1;
            ch1_concat_qrs = p_struct.HG2.QRS_annot_good_qual.ch1;
        elseif ~flag_empty(2) && (~sum(any(isnan(p_struct.HG1.ECG_good_qual.ch1))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch1))))
            ch1_concat_ecg = p_struct.HG1.ECG_good_qual.ch1;
            ch1_concat_qrs = p_struct.HG1.QRS_annot_good_qual.ch1;
        elseif flag_empty(2) && flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch1))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch1))))
            ch1_concat_ecg = NaN;
            ch1_concat_qrs = NaN;
        else
            ch1_concat_ecg = [p_struct.HG1.ECG_good_qual.ch1, p_struct.HG2.ECG_good_qual.ch1];
            ch1_concat_qrs = [p_struct.HG1.QRS_annot_good_qual.ch1, p_struct.HG2.QRS_annot_good_qual.ch1];
        end        

        p_struct.HG_concat.ECG_good_qual.ch1 = ch1_concat_ecg;
        p_struct.HG_concat.QRS_annot_good_qual.ch1 = ch1_concat_qrs;
    end
    if ~flag_empty(2) || ~flag_empty(3)
        if ~flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch2))) && ~sum(any(isnan(p_struct.HG2.ECG_good_qual.ch2))))
            ch2_concat_ecg = p_struct.HG2.ECG_good_qual.ch2;
            ch2_concat_qrs = p_struct.HG2.QRS_annot_good_qual.ch2;
        elseif ~flag_empty(2) && (~sum(any(isnan(p_struct.HG1.ECG_good_qual.ch2))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch2))))
            ch2_concat_ecg = p_struct.HG1.ECG_good_qual.ch2;
            ch2_concat_qrs = p_struct.HG1.QRS_annot_good_qual.ch2;
        elseif flag_empty(2) && flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch2))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch2))))
            ch2_concat_ecg = NaN;
            ch2_concat_qrs = NaN;
        else
            ch2_concat_ecg = [p_struct.HG1.ECG_good_qual.ch2, p_struct.HG2.ECG_good_qual.ch2];
            ch2_concat_qrs = [p_struct.HG1.QRS_annot_good_qual.ch2, p_struct.HG2.QRS_annot_good_qual.ch2];
        end

        p_struct.HG_concat.ECG_good_qual.ch2 = ch2_concat_ecg;
        p_struct.HG_concat.QRS_annot_good_qual.ch2 = ch2_concat_qrs;
    end

    if ~flag_empty(2) || ~flag_empty(3)   
        if ~flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch3))) && ~sum(any(isnan(p_struct.HG2.ECG_good_qual.ch3))))
            ch3_concat_ecg = p_struct.HG2.ECG_good_qual.ch3;
            ch3_concat_qrs = p_struct.HG2.QRS_annot_good_qual.ch3;
        elseif ~flag_empty(2) && (~sum(any(isnan(p_struct.HG1.ECG_good_qual.ch3))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch3))))
            ch3_concat_ecg = p_struct.HG1.ECG_good_qual.ch3;
            ch3_concat_qrs = p_struct.HG1.QRS_annot_good_qual.ch3;
        elseif flag_empty(2) && flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch3))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch3))))
            ch3_concat_ecg = NaN;
            ch3_concat_qrs = NaN;
        else
            ch3_concat_ecg = [p_struct.HG1.ECG_good_qual.ch3, p_struct.HG2.ECG_good_qual.ch3];
            ch3_concat_qrs = [p_struct.HG1.QRS_annot_good_qual.ch3, p_struct.HG2.QRS_annot_good_qual.ch3];
        end

        p_struct.HG_concat.ECG_good_qual.ch3 = ch3_concat_ecg;
        p_struct.HG_concat.QRS_annot_good_qual.ch3 = ch3_concat_qrs;
    end

    if empty_fun(p_struct)
        p_struct.message = 'participant file does not exist';
    end
    p_struct.has_AgCl = boolean(~flag_empty(1));
    p_struct.has_HG1 = boolean(~flag_empty(2));
    p_struct.has_HG2 = boolean(~flag_empty(3));

    if flag_empty(1)
        no_AgCl = [no_AgCl participants(p)];
    end
    if flag_empty(2)
        no_HG1 = [no_HG1 participants(p)];
    end
    if flag_empty(3)
        no_HG2 = [no_HG2 participants(p)];
    end

    if flag_empty(2) && flag_empty(3)
        no_HG = [no_HG participants(p)];
    end

    if sum(flag_empty) == 3
        no_all = [no_all participants(p)];
    end
    new_data_struct.(['p',char(participants(p))]) = p_struct;
end

new_data_struct.fs_ecg = data_struct.fs_ecg;
if save_files
    if ~use_filters
        save(fullfile(save_dir,"ECG_HG_manually_cleaned_quality_dataset_no_filters.mat"),"new_data_struct","no_HG1","no_HG2","no_HG","no_AgCl","no_all");
    else
        save(fullfile(save_dir,"ECG_HG_manually_cleaned_quality_dataset_MA.mat"),"new_data_struct","no_HG1","no_HG2","no_HG","no_AgCl","no_all");
        save(fullfile(save_dir,"metrics_manually_cleaned_deviation_from_noise.mat"),"metrics","no_HG1","no_HG2","no_HG","no_AgCl","no_all");
    end
end

