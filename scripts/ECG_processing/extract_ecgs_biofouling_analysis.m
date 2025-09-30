%% ECG Data Extraction for Biofouling Experiment Analysis
% Extracts and processes ECG data from 5 consecutive PPHG recordings in a
% single participant

close all
clear;
%% Setup environment
[scripts_dir, ~, ~] = fileparts(pwd);
[root_dir, ~, ~] = fileparts(scripts_dir);
setup_script = fullfile(root_dir,'utils','setup_environment.m');
run(setup_script);

%% Path to the dataset goes here
data_path = 'C:\Users\giann\OneDrive\Desktop\ECG HG paper\ecgs\';
data_struct = struct();

%% Path where results will be saved
save_dir = fullfile(root_dir,'data');
save_files = true;

%% Parameters
signals = ["1st","2nd","3rd","4th","5th"];
num_leads = 3;
filename = 'KMTPatient_00.kha';
hp_cutoff = 0.05; %BW HP filter cutoff (Hz) for detrending
hp_order = 4; %BW HP filter order 
lp_cutoff = 90; %BW LP filter cutoff (Hz) - 100 is the Nyquist frequency
lp_order = 4;
notch_freq1 = 50; %Central bandstop frequency (UAE power line frequency)
notch_freq2 = 80; %Power line 1st harmonic artifact 
notch_bw = 2; 

for s = 1:length(signals)
    
    ecg_path = fullfile(data_path,signals(s),filename);
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
    end
    
    fs_ecg = signal(1).samplefrequency;
    %Current analysis is for 3 leads
    L1 = signal(1).data;
    L2 = signal(2).data;
    L3 = signal(3).data;    
    shape = size(L1);
    L1 = reshape(L1,shape(1)*shape(2),1);
    L2 = reshape(L2,shape(1)*shape(2),1);
    L3 = reshape(L3,shape(1)*shape(2),1);
    ECG = double([L1, L2, L3]);
    
    %Remove flatline at end of each recording
    %Then isolate 10-second segments for each recording 10 seconds for visualizing
    if s == 1
        ECG = ECG(1:38000,:);
        ECG_seg = ECG(31000:31000+10*fs_ecg - 1,:);
        t_axis = ((1:38000) - 1)*(1/fs_ecg)';
        t_axis_seg = t_axis(31000:31000+10*fs_ecg - 1)';
    elseif s == 2
        ECG = ECG(1:36000,:);
        ECG_seg = ECG(14000:14000+10*fs_ecg - 1,:);
        t_axis = ((1:36000) - 1)*(1/fs_ecg)';
        t_axis_seg = t_axis(14000:14000+10*fs_ecg - 1)';
    elseif s == 3   
        ECG = ECG(1:46000,:);
        ECG_seg = ECG(10000:10000+10*fs_ecg - 1,:);
        t_axis = ((1:46000) - 1)*(1/fs_ecg)';
        t_axis_seg = t_axis(10000:10000+10*fs_ecg - 1)';
    elseif s == 4   
        ECG_seg = ECG(15000:15000+10*fs_ecg - 1,:);
        t_axis = ((1:size(ECG,1)) - 1)*(1/fs_ecg)';
        t_axis_seg = t_axis(15000:15000+10*fs_ecg - 1)';
    elseif s == 5
        ECG = ECG(1:34000,:);
        ECG_seg = ECG(28000:28000+10*fs_ecg - 1,:);
        t_axis = ((1:34000) - 1)*(1/fs_ecg)';
        t_axis_seg = t_axis(28000:28000+10*fs_ecg - 1)';
    end
    
    %Apply preprocessing for denoising to estimate SNR
    [blp, alp] = butter(lp_order, lp_cutoff / (fs_ecg/2), 'low');
    ECG_lp = filtfilt(blp, alp, ECG);
    Wn1 = [notch_freq1 - notch_bw, notch_freq1 + notch_bw] / (fs_ecg/2);
    Wn2 = [notch_freq2 - notch_bw, notch_freq2 + notch_bw] / (fs_ecg/2);
    [bnotch1, anotch1] = butter(2, Wn1, 'stop');
    [bnotch2, anotch2] = butter(2, Wn2, 'stop');
    ECG_notch1 = filtfilt(bnotch1, anotch1, ECG_lp);
    ECG_filtered = filtfilt(bnotch2, anotch2, ECG_notch1);
    %ECG_filtered = detrend_ecg_highpass(ECG_filtered,fs_ecg,hp_cutoff,hp_order);
    %ECG_filtered = detrend_ecg_movingavg(ECG_notch2,fs_ecg); %1-second window size

    noise_estimate = ECG - ECG_filtered;
    
    
    ECG_struct(s).ECG = ECG; 
    ECG_struct(s).ECG_seg = ECG_seg;
    ECG_struct(s).t_axis = t_axis;
    ECG_struct(s).t_axis_seg = t_axis_seg;

    %Calculate SNR
    SNR = zeros(3,1);
    for i = 1:3
        SNR(i)= snr(ECG(:,i), noise_estimate(:,i));
    end
    ECG_struct(s).SNR = SNR;
end

jsonStr = jsonencode(ECG_struct);
fid = fopen(fullfile(save_dir,'ECG_biofouling_data.json'),'w');
fwrite(fid,jsonStr,'char');
fclose(fid);
save(fullfile(save_dir,"ECG_biofouling_data.mat"),"ECG_struct")