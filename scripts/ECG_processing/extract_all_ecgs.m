%% Batch ECG Data Extraction
% Extracts and processes ECG data from all available participants
% Handles multiple electrode types (AgCl, HG) and applies quality filters
% Calculates signal quality metrics and saves processed data for further analysis

close all
clear;
%% Setup environment
[scripts_dir, ~, ~] = fileparts(pwd);
[root_dir, ~, ~] = fileparts(scripts_dir);
setup_script = fullfile(root_dir,'utils','setup_environment.m');
run(setup_script);

%% Path to the dataset goes here
data_path = 'C:\Users\giann\OneDrive\Desktop\ECG HG paper\participants_v2\';
data_struct = struct();

%% Path where results will be saved
save_dir = 'C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\';
save_files = true;

%% Parameters
signals = ["AgCl","HG1","HG2"];
num_leads = 3;
filename = 'KMTPatient_00.kha';
data_struct.fs_ecg = 200; %ecg sampling rate
data_struct.fs_imp = 1; %impedance sampling rate
window_size = data_struct.fs_ecg/data_struct.fs_imp;
num_windows = 2; %this defines the length of the segment
hp_cutoff = 0.05; %BW HP filter cutoff (Hz) for detrending
hp_order = 4; %BW HP filter order 
lp_cutoff = 90; %BW LP filter cutoff (Hz) - 100 is the Nyquist frequency
lp_order = 4;
notch_freq1 = 50; %Central bandstop frequency (UAE power line frequency)
notch_freq2 = 80; %Power line 1st harmonic artifact 
notch_bw = 2; 
thres_Imp = 100; %impedance > 100 kΩ considered not acceptable
use_filters = 1; %boolean, whether to use any filters in the analysis

%% 
participants = [];
for i = 1:43
    if i == 35 || i == 36 || i == 24
        continue
    end
    pid = string(i);
    participants = [participants,pid];
end

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
        
        ecg_path = fullfile(data_path,participants(p),signals(s),filename);
        if exist(ecg_path)
            try
                [header, signal,annotations]=ConvertKHA2MAT(ecg_path);
            catch ME
                if ((strcmp(ME.identifier,'MATLAB:badsubscript')) && length(ME.stack) == 4)
                    flag_empty(s) = 1;
                    for ch = 1:3
                        p_struct.(signal_label).ECG_good_qual.(['ch',int2str(ch)]) = NaN;
                        p_struct.(signal_label).Imp_good_qual.(['ch',int2str(ch)]) = NaN;
                    end
                    continue
                else
                    fprintf("Unknown error: %s \n",ME.identifier) 
                end
            end
        else
            flag_empty(s) = 1;
            for ch = 1:num_leads
                p_struct.(signal_label).ECG_good_qual.(['ch',int2str(ch)]) = NaN;
                p_struct.(signal_label).Imp_good_qual.(['ch',int2str(ch)]) = NaN;
            end
            continue
        end
        
        %Current analysis is for 3 leads
        L1 = signal(1).data;
        L2 = signal(2).data;
        L3 = signal(3).data;    

        Imp_L1 = signal(7).data;
        Imp_L2 = signal(8).data;
        Imp_L3 = signal(9).data;
        shape = size(L1);
        shape_imp = size(Imp_L1);
        L1 = reshape(L1,shape(1)*shape(2),1);
        L2 = reshape(L2,shape(1)*shape(2),1);
        L3 = reshape(L3,shape(1)*shape(2),1);
        ECG = double([L1, L2, L3]);
        
        Imp_L1 = double(reshape(Imp_L1,shape_imp(1)*shape_imp(2),1));
        avg_I1 = mean(abs(Imp_L1)); std_I1 = std(abs(Imp_L1));
        fprintf("Avg. |Imp| of Lead1: %d \n",avg_I1)
        fprintf("Std. |Imp| of Lead1: %d \n",std_I1)
        
        Imp_L2 = double(reshape(Imp_L2,shape_imp(1)*shape_imp(2),1));
        avg_I2 = mean(abs(Imp_L2)); std_I2 = std(abs(Imp_L2));
        fprintf("Avg. |Imp| of Lead2: %d \n",avg_I2)
        fprintf("Std. |Imp| of Lead2: %d \n",std_I2)

        Imp_L3 = double(reshape(Imp_L3,shape_imp(1)*shape_imp(2),1));
        avg_I3 = mean(abs(Imp_L3)); std_I3 = std(abs(Imp_L3));
        fprintf("Avg. |Imp| of Lead3: %d \n",avg_I3)
        fprintf("Std. |Imp| of Lead3: %d \n",std_I3)
        Imp = [Imp_L1, Imp_L2, Imp_L3];
          
        
        if use_filters
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

        %% Below is used to produce Figs 4c and S6 for participant p5 - p = 5
        if p == 5
            if s==1
                agcl = ECG(35000:41000,:);
                agcl_lp = ECG_lp(35000:41000,:);
                agcl_notch1 = ECG_just_notch1(35000:41000,:);
                agcl_notch2 = ECG_just_notch2(35000:41000,:);
                agcl_hp_VLF = ECG_just_hp_VLF(35000:41000,:); 
                agcl_MA = ECG_just_MA(35000:41000,:); 
                agcl_all_filters = ECG_filtered(35000:41000,:);  
                t_vec = linspace(0,size(agcl,1)/data_struct.fs_ecg,size(agcl,1));
                save(fullfile(root_dir,'data','agcl_p5.mat'),'t_vec','agcl','agcl_lp','agcl_notch1','agcl_notch2','agcl_hp_VLF','agcl_MA','agcl_all_filters')
            elseif s == 3
                hg = ECG(49000:55000,:);
                hg_lp = ECG_lp(49000:55000,:);
                hg_notch1 = ECG_just_notch1(49000:55000,:);
                hg_notch2 = ECG_just_notch2(49000:55000,:);
                hg_hp_VLF = ECG_just_hp_VLF(49000:55000,:); 
                hg_MA = ECG_just_MA(49000:55000,:); 
                hg_all_filters = ECG_filtered(49000:55000,:); 
                t_vec = linspace(0,size(hg,1)/data_struct.fs_ecg,size(hg,1));
                save(fullfile(root_dir,'data','hg_p5.mat'),'t_vec','hg','hg_lp','hg_notch1','hg_notch2','hg_hp_VLF','hg_MA','hg_all_filters')
            end
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
        if strcmp(participants(p),'1') || strcmp(participants(p),'5') || strcmp(participants(p),'8') || strcmp(participants(p),'9') || strcmp(participants(p),'10') || strcmp(participants(p),'13') || strcmp(participants(p),'23') || strcmp(participants(p),'25') || strcmp(participants(p),'41')
            plot_flag = true;
        end
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
            all_inds = find(abs(Imp(:,ch))<thres_Imp);
            if isempty(all_inds)
                ECG_good_qual.(['ch',int2str(ch)]) = NaN;
                imp_good_qual.(['ch',int2str(ch)]) = NaN;
                QRS_annot_good_qual.(['ch',int2str(ch)]) = NaN;
                fprintf("Discarding whole channel %d, of electrode %s, participant ID %s, because of high impedance... \n",ch,signals(s),participants(p))
                continue
            end
            inds = all_inds(1);
            %Ensure consecutive indices will not be kept
            for i = 2:length(all_inds)
                if abs(all_inds(i) - inds(end)) > 1
                    inds(end+1) = all_inds(i);  
                end
            end
            inds_start = max((inds-1)*(data_struct.fs_ecg/data_struct.fs_imp),1);
            inds_start(2:end) = inds_start(2:end) + 1;
            inds_end = inds_start + num_windows*window_size - 1;

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
            imp_good_qual.(['ch',int2str(ch)]) = Imp(inds(1:end-inds_to_remove));
            ECG_good_qual.(['ch',int2str(ch)]) = ecg_good_qual;
            QRS_annot_good_qual.(['ch',int2str(ch)]) = qrs_annot_good_qual(1:end-inds_to_remove);
            ECG_good_qual.(['ch',int2str(ch),'breakpoints']) = breakpoints;       
        end
        
        p_struct.(signal_label).ECG = ECG;
        p_struct.(signal_label).ECG_good_qual = ECG_good_qual;
        p_struct.(signal_label).QRS_annot_good_qual = QRS_annot_good_qual;
        p_struct.(signal_label).Imp = Imp;
        p_struct.(signal_label).Imp_good_qual = imp_good_qual;
        
    end
    if strcmp(participants(p),'1') || strcmp(participants(p),'5') || strcmp(participants(p),'8') || strcmp(participants(p),'9') || strcmp(participants(p),'10') || strcmp(participants(p),'13') || strcmp(participants(p),'23') || strcmp(participants(p),'25') || strcmp(participants(p),'41')
        plot_flag = true;
    else
        plot_flag = false;
    end
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
            ch1_concat_imp = p_struct.HG2.Imp_good_qual.ch1;
        elseif ~flag_empty(2) && (~sum(any(isnan(p_struct.HG1.ECG_good_qual.ch1))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch1))))
            ch1_concat_ecg = p_struct.HG1.ECG_good_qual.ch1;
            ch1_concat_qrs = p_struct.HG1.QRS_annot_good_qual.ch1;
            ch1_concat_imp = p_struct.HG1.Imp_good_qual.ch1;
        elseif flag_empty(2) && flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch1))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch1))))
            ch1_concat_ecg = NaN;
            ch1_concat_qrs = NaN;
            ch1_concat_imp = NaN;
        else
            ch1_concat_ecg = [p_struct.HG1.ECG_good_qual.ch1, p_struct.HG2.ECG_good_qual.ch1];
            ch1_concat_qrs = [p_struct.HG1.QRS_annot_good_qual.ch1, p_struct.HG2.QRS_annot_good_qual.ch1];
            ch1_concat_imp = [p_struct.HG1.Imp_good_qual.ch1'; p_struct.HG2.Imp_good_qual.ch1']';
        end        

        p_struct.HG_concat.ECG_good_qual.ch1 = ch1_concat_ecg;
        p_struct.HG_concat.QRS_annot_good_qual.ch1 = ch1_concat_qrs;
        p_struct.HG_concat.Imp_good_qual.ch1 = ch1_concat_imp;
    end
    if ~flag_empty(2) || ~flag_empty(3)
        if ~flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch2))) && ~sum(any(isnan(p_struct.HG2.ECG_good_qual.ch2))))
            ch2_concat_ecg = p_struct.HG2.ECG_good_qual.ch2;
            ch2_concat_qrs = p_struct.HG2.QRS_annot_good_qual.ch2;
            ch2_concat_imp = p_struct.HG2.Imp_good_qual.ch2;
        elseif ~flag_empty(2) && (~sum(any(isnan(p_struct.HG1.ECG_good_qual.ch2))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch2))))
            ch2_concat_ecg = p_struct.HG1.ECG_good_qual.ch2;
            ch2_concat_qrs = p_struct.HG1.QRS_annot_good_qual.ch2;
            ch2_concat_imp = p_struct.HG1.Imp_good_qual.ch2;
        elseif flag_empty(2) && flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch2))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch2))))
            ch2_concat_ecg = NaN;
            ch2_concat_qrs = NaN;
            ch2_concat_imp = NaN;
        else
            ch2_concat_ecg = [p_struct.HG1.ECG_good_qual.ch2, p_struct.HG2.ECG_good_qual.ch2];
            ch2_concat_qrs = [p_struct.HG1.QRS_annot_good_qual.ch2, p_struct.HG2.QRS_annot_good_qual.ch2];
            ch2_concat_imp = [p_struct.HG1.Imp_good_qual.ch2'; p_struct.HG2.Imp_good_qual.ch2']';
        end

        p_struct.HG_concat.ECG_good_qual.ch2 = ch2_concat_ecg;
        p_struct.HG_concat.QRS_annot_good_qual.ch2 = ch2_concat_qrs;
        p_struct.HG_concat.Imp_good_qual.ch2 = ch2_concat_imp;
    end
    
    if ~flag_empty(2) || ~flag_empty(3)   
        if ~flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch3))) && ~sum(any(isnan(p_struct.HG2.ECG_good_qual.ch3))))
            ch3_concat_ecg = p_struct.HG2.ECG_good_qual.ch3;
            ch3_concat_qrs = p_struct.HG2.QRS_annot_good_qual.ch3;
            ch3_concat_imp = p_struct.HG2.Imp_good_qual.ch3;
        elseif ~flag_empty(2) && (~sum(any(isnan(p_struct.HG1.ECG_good_qual.ch3))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch3))))
            ch3_concat_ecg = p_struct.HG1.ECG_good_qual.ch3;
            ch3_concat_qrs = p_struct.HG1.QRS_annot_good_qual.ch3;
            ch3_concat_imp = p_struct.HG1.Imp_good_qual.ch3;
        elseif flag_empty(2) && flag_empty(3) && (sum(any(isnan(p_struct.HG1.ECG_good_qual.ch3))) && sum(any(isnan(p_struct.HG2.ECG_good_qual.ch3))))
            ch3_concat_ecg = NaN;
            ch3_concat_qrs = NaN;
            ch3_concat_imp = NaN;
        else
            ch3_concat_ecg = [p_struct.HG1.ECG_good_qual.ch3, p_struct.HG2.ECG_good_qual.ch3];
            ch3_concat_qrs = [p_struct.HG1.QRS_annot_good_qual.ch3, p_struct.HG2.QRS_annot_good_qual.ch3];
            ch3_concat_imp = [p_struct.HG1.Imp_good_qual.ch3'; p_struct.HG2.Imp_good_qual.ch3']';
        end

        p_struct.HG_concat.ECG_good_qual.ch3 = ch3_concat_ecg;
        p_struct.HG_concat.QRS_annot_good_qual.ch3 = ch3_concat_qrs;
        p_struct.HG_concat.Imp_good_qual.ch3 = ch3_concat_imp;
    end
    
    if empty_fun(p_struct)
        p_struct.message = 'participant file does not exist';
    end
    p_struct.has_AgCl = boolean(~flag_empty(1));
    p_struct.has_HG1 = boolean(~flag_empty(2));
    p_struct.has_HG2 = boolean(~flag_empty(3));

    %% Log which participants didn't have which channels
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
    data_struct.(['p',char(participants(p))]) = p_struct;
end

if save_files
    if ~use_filters
        save(fullfile(save_dir,"ECG_HG_quality_dataset_no_filters.mat"),"data_struct","no_HG1","no_HG2","no_HG","no_AgCl","no_all");
    else
        save(fullfile(save_dir,"ECG_HG_quality_dataset_MA.mat"),"data_struct","no_HG1","no_HG2","no_HG","no_AgCl","no_all");
        save(fullfile(save_dir,"metrics_deviation_from_noise.mat"),"metrics","no_HG1","no_HG2","no_HG","no_AgCl","no_all");
    end
end

