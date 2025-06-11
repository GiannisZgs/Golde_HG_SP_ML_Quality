%We need to analyze all ECG channels as we are interested on the
%morphology - different waves are more pronounced across different leads
%Examine R peak similarity from ECG1 
%Examine T-wave similarity from ECG2 or 3
%Check if lead 1 gives consistently more prominent R peaks
%Same for leads 2,3
%Impedance can help us keep segments of good quality
%On the other hand, we need impedance to showcase the paper's argument
%regarding quality

%Participants not able to be accessed: 6(no files),11(cant load)
%P17 AgCl: all zero in the beginning, very high noise at second-half

clear;

dat_path = 'C:\Users\giann\OneDrive\Desktop\ECG HG paper\participants\';
data_struct = struct();
data_struct.fs_ecg = 200;
data_struct.fs_imp = 1;
relative_window_size = 2; %with respect to sampling periods e.g. 2 -> 2-sampling-periods-long window
window_size = data_struct.fs_ecg/data_struct.fs_imp;

thres_Imp = 100; %impedance > 100 kΩ considered not acceptable

participants = [];
for i = 2:42
    if i == 35 || i == 36
        continue
    end
    pid = string(i);
    participants = [participants,pid];
end

signals = ["Agcl","hydrogel 1","hydrogel 2"];

filename = 'KMTPatient_00.kha';

empty_fun = @(s) all(structfun(@isempty,s));%function to check if structure contains anything
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
        
        ecg_path = fullfile(dat_path,participants(p),signals(s),filename);
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
            for ch = 1:3
                p_struct.(signal_label).ECG_good_qual.(['ch',int2str(ch)]) = NaN;
                p_struct.(signal_label).Imp_good_qual.(['ch',int2str(ch)]) = NaN;
            end
            continue
        end
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
        
        %detrended_ECG_hp = detrend_ecg_highpass(ECG,data_struct.fs_ecg);
        detrended_ECG_MA = detrend_ecg_movingavg(ECG,data_struct.fs_ecg); %1-second window size
        
        %% Apply a QRS detection algorithm on Channel 1
        try
            [QRS_indexes,~] = pantompkins_qrs(detrended_ECG_MA(:,1),data_struct.fs_ecg);
        catch ME
            if strcmp(ME.identifier,'MATLAB:unassignedOutputs')
                QRS_indexes = [];
                fprintf("\n This recording will be flagged for removal - QRS detection failed because of no signal in part or whole of the recording \n")
            else
                rethrow(ME)
            end

        end
        %plot_pantompkins_qrs_detection(detrended_ECG_MA(:,1),filt_dat,int_dat,QRS_indexes,thF1,thI1);

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
            inds_end = inds_start + relative_window_size*(data_struct.fs_ecg/data_struct.fs_imp) - 1;

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
                ecg_good_qual = [ecg_good_qual detrended_ECG_MA(win_start:win_end,ch)];
                breakpoints = [breakpoints [win_start;win_end-1]];
                %Find qrs indexes that lie within the [inds_start,inds_end] interval    
                qrs_higher = find(QRS_indexes>=win_start);
                qrs_lower = find(QRS_indexes<=win_end);
                qrs_common = intersect(qrs_higher,qrs_lower);
                qrs_current = QRS_indexes(qrs_common);
                if any(qrs_current > relative_window_size*window_size)
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

        %ch1_breakpoints = [p_struct.HG1.ECG_good_qual.ch1breakpoints, p_struct.HG2.ECG_good_qual.ch1breakpoints];
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

        %ch2_breakpoints = [p_struct.HG1.ECG_good_qual.ch2breakpoints, p_struct.HG2.ECG_good_qual.ch2breakpoints];
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

        %ch3_breakpoints = [p_struct.HG1.ECG_good_qual.ch3breakpoints, p_struct.HG2.ECG_good_qual.ch3breakpoints];
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

save("C:\Users\giann\OneDrive\Desktop\ECG HG paper\ECG_HG_quality_dataset.mat","data_struct","no_HG1","no_HG2","no_HG","no_AgCl","no_all");

% to save that as mat file

%For every subject, load all signals (AgCl, Hg1, Hg2) and store them
%aligned
%Apply basic preprocessing - detrending

