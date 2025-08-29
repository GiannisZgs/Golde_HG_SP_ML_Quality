clear;
close all;
% Search windows around R-peak (in milliseconds)
% These define the regions where the function will look for P and T waves.
% Adjusting these can significantly impact detection accuracy.
p_search_window_ms = 250; % Max time before R-peak to look for P-wave
t_search_window_ms = 350; % Max time after R-peak to look for T-wave
calcMetrics = 0; %boolean to control metric calculation
%% Tune this up_percentile

up_percentile = 95;
bottom_percentile = 1;
use_filters = 1;
manually_cleaned = 0;
if manually_cleaned
    if use_filters
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\ECG_HG_manually_cleaned_quality_dataset_MA.mat");
    else
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\ECG_HG_manually_cleaned_quality_dataset_no_filters.mat");
    end
    data.data_struct = data.new_data_struct;
else
    if use_filters
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\ECG_HG_quality_dataset_MA.mat");
    else
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\ECG_HG_quality_dataset_no_filters.mat");
    end
end
try
    fs = data.data_struct.fs_ecg;
catch
    fprintf('Fs not found in data_struct \n')
    fs = 200;
end
p_search_window = (p_search_window_ms/1000)*fs;
t_search_window = (t_search_window_ms/1000)*fs;


no_AgCl = data.no_AgCl;
no_HG = data.no_HG;
no_AgCl_HG = [no_AgCl no_HG];
no_AgCl_HG_num = [];
for p = 1:length(no_AgCl_HG)
    no_AgCl_HG_num = [no_AgCl_HG_num str2num(no_AgCl_HG(p))];
end
no_AgCl_HG = unique(no_AgCl_HG_num);
data_struct = data.data_struct;
%Keep only participants that have: 1) AgCl, AND 2) At least one hydrogel
struct_fields = fields(data_struct);
for f = 1:length(struct_fields)
    fieldd = struct_fields(f); fieldd = fieldd{1};
    if strcmp(fieldd,"fs_ecg") || strcmp(fieldd,'fs_imp')
        continue
    end
    if any(no_AgCl_HG==str2num(fieldd(2:end)))
        fprintf("Skipping participant %s for incomplete data",fieldd)
        continue
    end
    
    p_struct = data_struct.(fieldd);
    AgCl = p_struct.AgCl;
    HG = p_struct.HG_concat;
    AgCl_ecg = AgCl.ECG_good_qual;
    AgCl_qrs_annotations = AgCl.QRS_annot_good_qual;
    HG_ecg = HG.ECG_good_qual;
    HG_qrs_annotations = HG.QRS_annot_good_qual;
            
    %Extract channels    
    
    agcl_segmented_beats = struct();
    for ch = 1:length(fields(AgCl_qrs_annotations))
        agcl_ecg = AgCl_ecg.(['ch',num2str(ch)]);
        agcl_annot = AgCl_qrs_annotations.(['ch',num2str(ch)]);  
        if any(isnan(agcl_ecg))
            agcl_segmented_beats.(['ch',num2str(ch)]) = NaN;
            continue
        end
        max_expected_heartbeats = 5*length(agcl_annot);
        heartbeats = cell(1,max_expected_heartbeats);
        c = 1; %count how many heartbeats we save
        for s = 1:size(agcl_ecg,2)
            agcl_ecg_segment = agcl_ecg(:,s);
            try
                agcl_qrs_segment = agcl_annot{1,s};
            catch
                agcl_qrs_segment = agcl_annot(1,s);
            end
            p_onsets = [];
            t_offsets = [];
            for r = 1:length(agcl_qrs_segment)
                rpeak = agcl_qrs_segment(r);
                if (rpeak-p_search_window <= 0) || (rpeak + t_search_window > length(agcl_ecg_segment)) 
                    p_onsets(end+1) = -1;
                    t_offsets(end+1) = -1;
                    agcl_qrs_segment(r) = -1;
                    continue 
                else
                    %[p_onset_idx, t_offset_idx] = detectPOnsetTOffset(hg_ecg_segment, fs, [rpeak],p_search_window_ms,t_search_window_ms);
                    %Use a simplified p_onset t_offset detection scheme for now
                    p_onsets(end+1) = rpeak - p_search_window;
                    t_offsets(end+1) = rpeak + t_search_window;
                end  
            end
            
            %Further segment the ECG wrt the found P-onsets & T-offsets
            rpeaks_to_keep_inds = find(agcl_qrs_segment > 0);
            if length(rpeaks_to_keep_inds) == 0
                continue
            end
            rpeaks = agcl_qrs_segment(rpeaks_to_keep_inds);
            p_onsets = p_onsets(rpeaks_to_keep_inds);
            t_offsets = t_offsets(rpeaks_to_keep_inds);
            
            for r = 1:length(rpeaks)
                heartbeats{c} = agcl_ecg_segment(p_onsets(r):t_offsets(r));
                c = c+1;
            end    
        end
        heartbeats = heartbeats(1:c-1);
        lengths = cellfun(@length,heartbeats);
        if length(unique(lengths)) == 1
            heartbeats = cell2mat(heartbeats);
        end
        if length(heartbeats) == 0
            agcl_segmented_beats.(['ch',num2str(ch)]) = NaN;
            continue
        end
        %% Outlier detection
        %Exclude top5 and bottom5 percentiles
        variance_heartbeats = var(heartbeats);
        var_sorted = sort(variance_heartbeats);
        lower_thres = prctile(var_sorted,bottom_percentile);
        upper_thres = prctile(var_sorted,up_percentile);
        inds_to_keep = (variance_heartbeats >= lower_thres) & (variance_heartbeats <= upper_thres);
        heartbeats = heartbeats(:,inds_to_keep);
        agcl_segmented_beats.(['ch',num2str(ch)]) = heartbeats;
    end      
    
%     if ~any(isnan(agcl_segmented_beats.ch1))
%         figure;
%         plot(agcl_segmented_beats.ch1)
%         title('AgCl Channel 1 Heartbeat Profile')
%     end
%     
%     if ~any(isnan(agcl_segmented_beats.ch2))
%         figure;
%         plot(agcl_segmented_beats.ch2)
%         title('AgCl Channel 2 Heartbeat Profile')
%     end
%     
%     if ~any(isnan(agcl_segmented_beats.ch3))
%         figure;
%         plot(agcl_segmented_beats.ch3)
%         title('AgCl Channel 3 Heartbeat Profile')
%     end    
    
    %% Hydrogel
    hg_segmented_beats = struct();
    for ch = 1:length(fields(HG_qrs_annotations))
        
        hg_ecg = HG_ecg.(['ch',num2str(ch)]);
        hg_annot = HG_qrs_annotations.(['ch',num2str(ch)]);
        if any(isnan(hg_ecg))
            hg_segmented_beats.(['ch',num2str(ch)]) = NaN;
            continue
        end
        max_expected_heartbeats = 5*length(hg_annot);
        heartbeats = cell(1,max_expected_heartbeats);
        c = 1; %count how many heartbeats we save
        for s = 1:size(hg_ecg,2)
            hg_ecg_segment = hg_ecg(:,s);
            try
                hg_qrs_segment = hg_annot{1,s};
            catch
                hg_qrs_segment = hg_annot(1,s);
            end
            %rri_hg = diff(hg_qrs_segment);
            p_onsets = [];
            t_offsets = [];
            for r = 1:length(hg_qrs_segment)
                rpeak = hg_qrs_segment(r);
                if (rpeak-p_search_window <= 0) || (rpeak + t_search_window > length(hg_ecg_segment)) 
                    p_onsets(end+1) = -1;
                    t_offsets(end+1) = -1;
                    hg_qrs_segment(r) = -1;
                    continue 
                else
                    %[p_onset_idx, t_offset_idx] = detectPOnsetTOffset(hg_ecg_segment, fs, [rpeak],p_search_window_ms,t_search_window_ms);
                    %Use a simplified p_onset t_offset detection scheme for now
                    p_onsets(end+1) = rpeak - p_search_window;
                    t_offsets(end+1) = rpeak + t_search_window;
                end  
            end
            
            %Further segment the ECG wrt the found P-onsets & T-offsets
            rpeaks_to_keep_inds = find(hg_qrs_segment > 0);
            if length(rpeaks_to_keep_inds) == 0
                continue
            end
            rpeaks = hg_qrs_segment(rpeaks_to_keep_inds);
            p_onsets = p_onsets(rpeaks_to_keep_inds);
            t_offsets = t_offsets(rpeaks_to_keep_inds);
            
            for r = 1:length(rpeaks)
                heartbeats{c} = hg_ecg_segment(p_onsets(r):t_offsets(r));
                c = c+1;
            end    
        end
        heartbeats = heartbeats(1:c-1);
        lengths = cellfun(@length,heartbeats);
        if length(unique(lengths)) == 1
            heartbeats = cell2mat(heartbeats);
        end
        if length(heartbeats) == 0
            hg_segmented_beats.(['ch',num2str(ch)]) = NaN;
            continue
        end
        %% Outlier detection
        %Exclude top5 and bottom5 percentiles
        variance_heartbeats = var(heartbeats);
        var_sorted = sort(variance_heartbeats);
        lower_thres = prctile(var_sorted,bottom_percentile);
        upper_thres = prctile(var_sorted,up_percentile);
        inds_to_keep = (variance_heartbeats >= lower_thres) & (variance_heartbeats <= upper_thres);
        heartbeats = heartbeats(:,inds_to_keep);
        hg_segmented_beats.(['ch',num2str(ch)]) = heartbeats;
        
    end
    
    if ~any(isnan(hg_segmented_beats.ch1))
        figure;
        plot(hg_segmented_beats.ch1)
        title('HG Channel 1 Heartbeat Profile')
    end
    
    if ~any(isnan(hg_segmented_beats.ch2))
        figure;
        plot(hg_segmented_beats.ch2)
        title('HG Channel 2 Heartbeat Profile')
    end
    
    
    if ~any(isnan(hg_segmented_beats.ch3))
        figure;
        plot(hg_segmented_beats.ch3)
        title('HG Channel 3 Heartbeat Profile')
    end
    
    %plot(agcl_ecg_segment)
    %hold on
    %plot(agcl_qrs_segment,agcl_ecg_segment(agcl_qrs_segment),'ok');
    %plot(hg_ecg_segment)
    %plot(hg_qrs_segment,hg_ecg_segment(hg_qrs_segment),'og');
    
    %% Save heartbeat profiles for post-processing
    profiling_struct.(fieldd).HG = hg_segmented_beats;
    profiling_struct.(fieldd).AgCl = agcl_segmented_beats;
    if calcMetrics
        %% Calculate quality metrics
        %Inter-sensor for same channel of different sensors (same placement)
        %Intra-sensor: differences between peaks of the same recording
        for ch = 1:length(fields(HG_qrs_annotations))

            %Inter-sensor Cross-correlation 
            if (~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)]))))) && ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)]))))))    
                [CorrMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]),LagMat.(struct_fields{f}).(['ch',num2str(ch)])] = batch_inter_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), agcl_segmented_beats.(['ch',num2str(ch)]), 'xcorr');
            else
                CorrMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]) = [];
                LagMat.(struct_fields{f}).(['ch',num2str(ch)]) = [];
            end
            %Intra-sensor correlation
            if ~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)])))))
                [corrmat,~] = batch_intra_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]),'xcorr');
                CorrMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = triu(corrmat,1); %Keep upper_triangular matrix
            else
                CorrMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
            end
            if ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))    
                [corrmat,~] = batch_intra_sensor_metrics(agcl_segmented_beats.(['ch',num2str(ch)]),'xcorr');
                CorrMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = triu(corrmat,1);
            else
                CorrMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
            end

            %Inter-sensor NRMSE
            if (~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)]))))) && ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))) 
                NrmseMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)])= batch_inter_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), agcl_segmented_beats.(['ch',num2str(ch)]), 'nrmse');
            else
                NrmseMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]) = [];
            end
            %Intra-sensor NRMSE
            if ~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)])))))
                nrmsematHG= batch_intra_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), 'nrmse');
                NrmseMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = triu(nrmsematHG,1);
            else
                NrmseMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
            end
            if ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))    
                nrmsematAgCl= batch_intra_sensor_metrics(agcl_segmented_beats.(['ch',num2str(ch)]), 'nrmse');
                NrmseMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = triu(nrmsematAgCl,1);
            else
                NrmseMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
            end

            %Inter-sensor JSD of Extreme Values
            if (~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)]))))) && ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)]))))))    
                [JsdMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]),~,p1_jsd_x,p2_jsd_x,edges_jsd_x] = batch_inter_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), agcl_segmented_beats.(['ch',num2str(ch)]), 'jsd',100);
            else
                JsdMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]) = [];
            end

            %Intra-sensor JSD of Extreme Values
            if ~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)])))))
                [jsdmatHG,~,p1_jsd_hg,p2_jsd_hg,edges_jsd_hg] = batch_intra_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), 'jsd',100);
                JsdMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = triu(jsdmatHG,1);
            else
                JsdMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
            end
            if ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))    
                [jsdmatAgCl,~,p1_jsd_ag,p2_jsd_ag,edges_jsd_ag] = batch_intra_sensor_metrics(agcl_segmented_beats.(['ch',num2str(ch)]), 'jsd',100);
                JsdMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = triu(jsdmatAgCl,1);
            else
                JsdMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
            end

            %Inter-sensor Cosine Similarity
            if (~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)]))))) && ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))) 
                CosMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]) = batch_inter_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), agcl_segmented_beats.(['ch',num2str(ch)]), 'cosine');    
            else
                CosMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]) = [];
            end

            %Intra-sensor Cosine Similarity
            if ~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)])))))
                cosmatHG = batch_intra_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), 'cosine');
                CosMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = triu(cosmatHG);
            else
                CosMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
            end
            if ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))    
                cosmatAgCl = batch_intra_sensor_metrics(agcl_segmented_beats.(['ch',num2str(ch)]), 'cosine');
                CosMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = triu(cosmatAgCl);
            else
                CosMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
            end

            %Inter-sensor Earth Mover's Distance
            if (~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)]))))) && ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)]))))))   
                [EmdMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]),p1_emd_x,p2_emd_x,edges_emd_x] = batch_inter_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), agcl_segmented_beats.(['ch',num2str(ch)]), 'emd',100);        
            else
                EmdMat.(struct_fields{f}).Xsensor.(['ch',num2str(ch)]) = [];
            end
            %Intra-sensor Earth Mover's Distance
            if ~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)])))))
                [emdmatHG,p1_emd_hg,p2_emd_hg] = batch_intra_sensor_metrics(hg_segmented_beats.(['ch',num2str(ch)]), 'emd',100);                    
                EmdMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = triu(emdmatHG,1);
            else
                EmdMat.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
            end
            if ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))    
                [emdmatAgCl,p1_emd_ag,p2_emd_ag] = batch_intra_sensor_metrics(agcl_segmented_beats.(['ch',num2str(ch)]), 'emd',100);
                EmdMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = triu(emdmatAgCl,1);
            else
                EmdMat.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
            end

            % Calculation of Signal quality (Only intra)
            if ~(any(any(isnan(agcl_segmented_beats.(['ch',num2str(ch)])))))    
                [H_w,H_s,ZCR1,ZCR2,SF,Compr] = waveform_quality_scores(agcl_segmented_beats.(['ch',num2str(ch)]));
                SigQual.WaveEntropy.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = H_w;
                SigQual.SpecEntropy.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = H_s;
                SigQual.ZCR1.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = ZCR1;
                SigQual.ZCR2.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = ZCR2;
                SigQual.SpecFlatness.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = SF;
                SigQual.Compress.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = Compr;
            else
                SigQual.WaveEntropy.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
                SigQual.SpecEntropy.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
                SigQual.ZCR1.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
                SigQual.ZCR2.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
                SigQual.SpecFlatness.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
                SigQual.Compress.(struct_fields{f}).AgCl.(['ch',num2str(ch)]) = [];
            end
            if ~(any(any(isnan(hg_segmented_beats.(['ch',num2str(ch)])))))
                [H_w,H_s,ZCR1,ZCR2,SF,Compr] = waveform_quality_scores(hg_segmented_beats.(['ch',num2str(ch)]));
                SigQual.WaveEntropy.(struct_fields{f}).HG.(['ch',num2str(ch)]) = H_w;
                SigQual.SpecEntropy.(struct_fields{f}).HG.(['ch',num2str(ch)]) = H_s;
                SigQual.ZCR1.(struct_fields{f}).HG.(['ch',num2str(ch)]) = ZCR1;
                SigQual.ZCR2.(struct_fields{f}).HG.(['ch',num2str(ch)]) = ZCR2;
                SigQual.SpecFlatness.(struct_fields{f}).HG.(['ch',num2str(ch)]) = SF;
                SigQual.Compress.(struct_fields{f}).HG.(['ch',num2str(ch)]) = Compr;
            else
                SigQual.WaveEntropy.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
                SigQual.SpecEntropy.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
                SigQual.ZCR1.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
                SigQual.ZCR2.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
                SigQual.SpecFlatness.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
                SigQual.Compress.(struct_fields{f}).HG.(['ch',num2str(ch)]) = [];
            end
        end
        %avg_corr_ch1 = mean(reshape(CorrMat_ch1,1,size(CorrMat_ch1,1)*size(CorrMat_ch1,2)));
        %std_corr_ch1 = std(reshape(CorrMat_ch1,1,size(CorrMat_ch1,1)*size(CorrMat_ch1,2)));
        %avg_corr_ch1_hg = mean(CorrMat_ch1_HG(find(CorrMat_ch1_HG~=0)));
        %std_corr_ch1_hg = std(CorrMat_ch1_HG(find(CorrMat_ch1_HG~=0)));    
        %avg_mse_ch1 = mean(reshape(MSE_mat_ch1,1,size(MSE_mat_ch1,1)*size(MSE_mat_ch1,2)));
        %std_mse_ch1 = std(reshape(MSE_mat_ch1,1,size(MSE_mat_ch1,1)*size(MSE_mat_ch1,2)));
        %avg_nrmse_ch1 = mean(reshape(NRMSE_mat_ch1,1,size(NRMSE_mat_ch1,1)*size(NRMSE_mat_ch1,2)));
        %std_nrmse_ch1 = std(reshape(NRMSE_mat_ch1,1,size(NRMSE_mat_ch1,1)*size(NRMSE_mat_ch1,2)));              
        %avg_cos_ch1 = mean(reshape(Cos_mat_ch1,1,size(Cos_mat_ch1,1)*size(Cos_mat_ch1,2)));
        %std_cos_ch1 = std(reshape(Cos_mat_ch1,1,size(Cos_mat_ch1,1)*size(Cos_mat_ch1,2)));
    end
    fprintf("Finished results extraction of subject %s \n",struct_fields{f}) 
end
%save('C:\Users\giann\OneDrive\Desktop\ECG HG paper\similarity_analysis_results.mat','CorrMat','NrmseMat','CosMat','JsdMat','EmdMat','SigQual','-v7.3');

%% Save profiles for further analysis in MATLAB
if manually_cleaned
    if use_filters
        save('C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\manually_cleaned_heartbeat_profiles_MA.mat','profiling_struct')
    else
        save('C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\manually_cleaned_heartbeat_profiles_no_filters.mat','profiling_struct')
    end
else 
    if use_filters
        save('C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\heartbeat_profiles_MA.mat','profiling_struct')
    else
        save('C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\heartbeat_profiles_no_filters.mat','profiling_struct')
    end
end
%% Save results for processing in R
if calcMetrics
    results_struct.CorrMat = CorrMat;
    results_struct.NrmseMat = NrmseMat;
    results_struct.CosMat = CosMat;
    results_struct.JsdMat = JsdMat;
    results_struct.EmdMat = EmdMat;
    results_struct.SigQual = SigQual;

    jsonStr = jsonencode(results_struct);
    fid = fopen('C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\similarity_analysis_results.json','w');
    fwrite(fid,jsonStr,'char');
    fclose(fid);
end