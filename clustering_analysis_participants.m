clear;
close all;
%% Implement a clustering analysis to extract average heartbeat profiles
%% This should be done across patients, across sensors, across channels
use_filters = 0;
plot_profiles = true;
manually_cleaned = 1;
if manually_cleaned
    if use_filters
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\manually_cleaned_heartbeat_profiles_MA.mat");
    else
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\manually_cleaned_heartbeat_profiles_no_filters.mat");
    end
else    
    if use_filters
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\heartbeat_profiles_MA.mat");
    else
        data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\heartbeat_profiles_no_filters.mat");
    end
end
profiles = data.profiling_struct;
num_clusters_multiplier = 1; %number of clusters will be the number of participants times this
leads = 3;
sensors = ["AgCl","HG"];
sensors_num = 2;
fs = 200;
rng(1); % for reproducibility
%% Take all channels of all participants, both sensors and cluster them into N_participants clusters

for s = 1:sensors_num
    for ch=1:leads

        channel = [];
        labels = [];
        participants = fields(profiles);
        participants_used = 0;
        for p = 1:length(participants)
            participant = participants(p); participant = participant{1};
            if ~any(any(isnan(profiles.(participant).(sensors(s)).(['ch',num2str(ch)]))))
                signal = profiles.(participant).(sensors(s)).(['ch',num2str(ch)]);
                channel = [channel signal];
                labels = [labels str2double(participant(2:end))*ones(1,size(signal,2))];
                participants_used = participants_used + 1;
                %centroids.(sensors(s)).(['ch',num2str(ch)]).(participant).direct_centroid_waveform = mean(signal,2);
            else
                fprintf('%s, %s, Ch. %d not used \n', participant, sensors(s), ch);
            end
        end

        channel = channel'; %(samples, features) arrangement
        
        %% Reduce data to visualize - PCA
        %Ensure kept components explain at least 99% of the variability
        %Don't go over 50 components
        [coeff, score, latent, ~, explained,mu] = pca(channel);
        cumulative_explained = cumsum(explained);
        num_components = find(cumulative_explained >= 99, 1, 'first');
        num_components = min(num_components, 50);
        reduced_data_99 = score(:,1:num_components);

        num_clusters = participants_used*num_clusters_multiplier;
        %[idx_kmeans_pca,centroids_kmeans_pca] = kmeans(reduced_data_99,num_clusters,'MaxIter',1000); %'Replicates',5,'Start','cluster','MaxIter',1000,'Distance','sqeuclidean');
        [idx_kmeans_pca,centroids_kmeans_pca] = kmeans(reduced_data_99,num_clusters,'Replicates',5,'Start','cluster','MaxIter',1000,'Distance','sqeuclidean');
        reduced_data_with_centroids = [reduced_data_99; centroids_kmeans_pca];

        %centroids_pca = (centroids_kmeans - mu) * coeff(:,1:num_components);
        %reduced_data_with_centroids = [reduced_data_99; centroids_pca];

        %% Reduce further to 2 dimensions for visualizing
        combined_tsne = tsne(reduced_data_with_centroids, 'NumDimensions', 2, 'Perplexity', 30);

        tsne_data = combined_tsne(1:end-num_clusters, :);
        tsne_centroids = combined_tsne(end-num_clusters+1:end, :);

        %% Inverse transform of centroids from principal component space to time domain space
        centroids_waveforms = centroids_kmeans_pca * coeff(:,1:num_components)' + mu;

        %% Compute clustering accuracy and map clusters to true labels - Assignments(j) is the true label corresponding to cluster j
        labels_encoded = grp2idx(labels); % convert to categorical integers
        unique_encoded = unique(labels_encoded);
        unique_labels = unique(labels);
        [acc_kmeans,macro_f1_kmeans,micro_f1_kmeans,cluster_assignments] = greedy_cluster_mapping(labels_encoded', idx_kmeans_pca);
        for j = 1:length(cluster_assignments)
            participant = ['p',num2str(unique_labels(cluster_assignments(j)))];% participant = participant{1};
            centroids.(sensors(s)).(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform = centroids_waveforms(j,:)';
        end
        fprintf('%s Ch. %d Unsupervised Participant Identification Accuracy (K-Means): %.2f%%\n', sensors(s), ch, acc_kmeans * 100);
        fprintf('%s Ch. %d Unsupervised Participant Identification Macro-F1 (K-Means): %.2f%%\n', sensors(s), ch, macro_f1_kmeans * 100);
        fprintf('%s Ch. %d Unsupervised Participant Identification Micro-F1 (K-Means): %.2f%%\n', sensors(s), ch, micro_f1_kmeans * 100);

%         participant_identification.(sensors(s)).(['ch',num2str(ch)]).acc = acc_kmeans;
%         participant_identification.(sensors(s)).(['ch',num2str(ch)]).macro_f1 = macro_f1_kmeans;
%         participant_identification.(sensors(s)).(['ch',num2str(ch)]).micro_f1 = micro_f1_kmeans;

        %% Visualize with colored clusters

        figure;
        subplot(1,2,1);
        gscatter(tsne_data(:,1), tsne_data(:,2), labels);
        title('Ground Truth - TSNE manifold');
        xlabel('t-SNE 1'); ylabel('t-SNE 2');

        subplot(1,2,2);
        gscatter(tsne_data(:,1), tsne_data(:,2), idx_kmeans_pca);
        hold on; scatter(tsne_centroids(:,1), tsne_centroids(:,2), 100, 'kx', 'LineWidth', 2);
        title('K-Means Clustering on PCA projection - TSNE manifold');
        xlabel('t-SNE 1'); ylabel('t-SNE 2');
        lgd = findobj('type','legend');
        delete(lgd)
        %legend(['Accuracy: ' sprintf('%.2f%%', acc_kmeans * 100)], 'Location', 'northeast');
        %suptitle = [char(sensors(s)) ' - Ch.' num2str(ch)];
        suptitle = [char(sensors(s)) ' - Ch.' num2str(ch) ' - Participant Identification - Macro F1: ' sprintf('%.2f%%', macro_f1_kmeans * 100) ' / Accuracy: ' sprintf('%.2f%%', acc_kmeans * 100)];
        sgtitle(suptitle)
        
        %% Store TSNE data used in plots for further processing
        tsne_vis.(sensors(s)).(['ch',num2str(ch)]).projection = tsne_data;
        tsne_vis.(sensors(s)).(['ch',num2str(ch)]).centroids = tsne_centroids;
        tsne_vis.(sensors(s)).(['ch',num2str(ch)]).kmeans_clusters = idx_kmeans_pca;
        tsne_vis.(sensors(s)).(['ch',num2str(ch)]).ground_truth = labels;
        tsne_vis.(sensors(s)).(['ch',num2str(ch)]).acc = acc_kmeans;
        tsne_vis.(sensors(s)).(['ch',num2str(ch)]).macro_f1 = macro_f1_kmeans;
        tsne_vis.(sensors(s)).(['ch',num2str(ch)]).micro_f1 = micro_f1_kmeans;
    
    end
end

%% Visualize some profiles for different subjects
%Cluster profiles

plot_count = 1;
if manually_cleaned
    plot_participants = {'p1','p5','p10','p39'};
else
    plot_participants = {'p1','p5','p8', 'p9','p10','p13','p23','p25','p41'};
    %plot_participants = fields(profiles);
end
    %figure;
for p = 1:length(plot_participants)
    participant = plot_participants(p); participant = participant{1};
    plot_count = 1;
    if plot_profiles
        figure;
    end
    global_min = 10000000;
    global_max = -10000000;
    for s = 1:2
        for ch = 1:3
            if ~(any(any(isnan(profiles.(participant).(sensors(s)).(['ch',num2str(ch)])))))
                cluster_profile = centroids.(sensors(s)).(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform;
                %subplot(length(plot_participants),6,plot_count)
                if plot_profiles
                    subplot(1,6,plot_count)
                    plot(cluster_profile,'k','LineWidth',1.5)
                    xlabel((['ch',num2str(ch)])); %
                    if ch==1 && s==1
                        ylabel(participant);
                    end
                    if ch == 2
                        title(sensors(s))
                    end
                    if min(cluster_profile) < global_min
                        global_min = min(cluster_profile);
                    end
                    if max(cluster_profile) > global_max
                        global_max = max(cluster_profile);
                    end
                end
                %subplot(length(plot_participants),6,plot_count)
                %subplot(1,6,plot_count)
                features.(participant).(sensors(s)).(['ch',num2str(ch)]) = extract_ecg_heartbeat_features(cluster_profile, fs);
            end
            plot_count = plot_count + 1;
            
        end
    end
    if plot_profiles
        for i = 1:6
            subplot(1,6,i)
            ylim([global_min,global_max])
        end
    end
end
%sgtitle('Centroid Profiles obtained from Clustering Analysis for 4 participants')

%Waveform profiles

% plot_count = 1;
% plot_participants = {'p2','p4','p5','p8'};
% figure;
% for p = 1:length(plot_participants)
%     participant = plot_participants(p); participant = participant{1};
%     for s = 1:2
%         for ch = 1:3
%             if ~(any(any(isnan(profiles.(participant).(sensors(s)).(['ch',num2str(ch)])))))
%                 waveform_avg_profile = centroids.(sensors(s)).(['ch',num2str(ch)]).(participant).direct_centroid_waveform;
%                 subplot(length(plot_participants),6,plot_count)
%                 plot(cluster_profile,'k','LineWidth',1.5)
%                 xlabel((['ch',num2str(ch)])); %
%                 if ch==1 && s==1
%                     ylabel(participant);
%                 end
%                 if ch == 2
%                     title(sensors(s))
%                 end
%                 subplot(length(plot_participants),6,plot_count)
%                 plot_count = plot_count + 1;
%             end
%         end
%     end
% end
% sgtitle('Average Waveform Profiles for 4 participants')
                

%% Calculate inter-sensor metrics for every channel
% xcorr11_from_cluster = [];
% mse11_from_cluster = [];
% xcorr22_from_cluster = [];
% mse22_from_cluster = [];
% xcorr33_from_cluster = [];
% mse33_from_cluster = [];
% xcorr11_simple_avg = [];
% mse11_simple_avg = [];
% xcorr22_simple_avg = [];
% mse22_simple_avg = [];
% xcorr33_simple_avg = [];
% mse33_simple_avg = [];
for p = 1:length(participants)
    participant = participants(p); participant = participant{1};
    for ch = 1:3
        if ~(any(any(isnan(profiles.(participant).AgCl.(['ch',num2str(ch)])))) || any(any(isnan(profiles.(participant).HG.(['ch',num2str(ch)])))))
            xcorr_cluster = batch_inter_sensor_metrics(centroids.AgCl.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,centroids.HG.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,'xcorr');
            %xcorr_simple_avg = batch_inter_sensor_metrics(centroids.AgCl.(['ch',num2str(ch)]).(participant).direct_centroid_waveform,centroids.HG.(['ch',num2str(ch)]).(participant).direct_centroid_waveform,'xcorr');
            mse_cluster = batch_inter_sensor_metrics(centroids.AgCl.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,centroids.HG.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,'nrmse');
            %mse_simple_avg = batch_inter_sensor_metrics(centroids.AgCl.(['ch',num2str(ch)]).(participant).direct_centroid_waveform,centroids.HG.(['ch',num2str(ch)]).(participant).direct_centroid_waveform,'nrmse');

            %Save in hierarchical structure
            xcorr.(['ch',num2str(ch)]).(participant) = xcorr_cluster;
            %xcorr.(['ch',num2str(ch)]).(participant).simple_avg = xcorr_simple_avg;
            mse.(['ch',num2str(ch)]).(participant) = mse_cluster;
            %mse.(['ch',num2str(ch)]).(participant).simple_avg = mse_simple_avg;
%             if ch == 1
%                 xcorr11_from_cluster = [xcorr11_from_cluster; xcorr_cluster];
%                 %xcorr11_simple_avg = [xcorr11_simple_avg; xcorr_simple_avg];
%                 mse11_from_cluster = [mse11_from_cluster; mse_cluster];
%                 %mse11_simple_avg = [mse11_simple_avg; mse_simple_avg];
%             elseif ch == 2
%                 xcorr22_from_cluster = [xcorr22_from_cluster; xcorr_cluster];
%                 %xcorr22_simple_avg = [xcorr22_simple_avg; xcorr_simple_avg];
%                 mse22_from_cluster = [mse22_from_cluster; mse_cluster];
%                 %mse22_simple_avg = [mse22_simple_avg; mse_simple_avg];
%             elseif ch == 3
%                 xcorr33_from_cluster = [xcorr33_from_cluster; xcorr_cluster];
%                 %xcorr33_simple_avg = [xcorr33_simple_avg; xcorr_simple_avg];
%                 mse33_from_cluster = [mse33_from_cluster; mse_cluster];
%                 %mse33_simple_avg = [mse33_simple_avg; mse_simple_avg];
%             end
        end
    end
end

%% Cross-correlation boxplots
% figure;
% subplot(1,3,1); %subplot(2,3,1);
% boxplot(xcorr11_from_cluster)
% ylim([0.2,1.1])
% %title('Clustering Centroid profiles AgCl-HG cross-correlation');
% xlabel('Channel 1'); ylabel('Cross-correlation');
% 
% subplot(1,3,2); %subplot(2,3,2);
% boxplot(xcorr22_from_cluster)
% ylim([0.2,1.1])
% title('Clustering Centroid profiles AgCl-HG cross-correlation');
% xlabel('Channel 2'); %ylabel('Cross-correlation');
% 
% subplot(1,3,3); %subplot(2,3,3);
% boxplot(xcorr33_from_cluster)
% ylim([0.2,1.1])
% %title('Clustering Centroid profiles AgCl-HG cross-correlation');
% xlabel('Channel 3'); %ylabel('Cross-correlation');
% 
% subplot(2,3,4);
% boxplot(xcorr11_simple_avg)
% ylim([0.2,1.1])
% %title('Simple Average profiles AgCl-HG cross-correlation');
% xlabel('Channel 1'); ylabel('Cross-correlation');
% 
% subplot(2,3,5);
% boxplot(xcorr22_simple_avg)
% ylim([0.2,1.1])
% title('Simple Average profiles AgCl-HG cross-correlation');
% xlabel('Channel 2'); %ylabel('Cross-correlation');
% 
% subplot(2,3,6);
% boxplot(xcorr33_simple_avg)
% ylim([0.2,1.1])
% %title('Simple Average profiles AgCl-HG cross-correlation');
% xlabel('Channel 3'); %ylabel('Cross-correlation');

%% Mean Squared Error boxplots
% figure;
% subplot(1,3,1); %subplot(1,3,1);
% boxplot(mse11_from_cluster)
% ylim([0,1])
% %title('Clustering Centroid profiles AgCl-HG cross-correlation');
% xlabel('Channel 1'); ylabel('Normalized Root Mean Squared Error');
% 
% subplot(1,3,2); %subplot(1,3,2);
% boxplot(mse22_from_cluster)
% ylim([0,1])
% title('Clustering Centroid profiles AgCl-HG Normalized Root Mean Squared Error');
% xlabel('Channel 2'); %ylabel('Cross-correlation');
% 
% subplot(1,3,3); %subplot(2,3,3);
% boxplot(mse33_from_cluster)
% ylim([0,1])
% %title('Clustering Centroid profiles AgCl-HG cross-correlation');
% xlabel('Channel 3'); %ylabel('Cross-correlation');
% 
% subplot(2,3,4);
% boxplot(mse11_simple_avg)
% ylim([0,1])
% %title('Simple Average profiles AgCl-HG cross-correlation');
% xlabel('Channel 1'); ylabel('Normalized Root Mean Squared Error');
% 
% subplot(2,3,5);
% boxplot(mse22_simple_avg)
% ylim([0,1])
% title('Simple Average profiles AgCl-HG Normalized Root Mean Squared Error');
% xlabel('Channel 2'); %ylabel('Cross-correlation');
% 
% subplot(2,3,6);
% boxplot(mse33_simple_avg)
% ylim([0,1])
% %title('Simple Average profiles AgCl-HG cross-correlation');
% xlabel('Channel 3'); %ylabel('Cross-correlation');

%% Save results and data for further processing
if manually_cleaned
    if use_filters
    save("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\manually_cleaned_participant_id_results.mat","centroids","features","tsne_vis","xcorr","mse");
else
    save("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\manually_cleaned_participant_id_results_no_filters.mat","centroids","features","tsne_vis","xcorr","mse");
    end
else
    if use_filters
        save("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\participant_id_results.mat","centroids","features","tsne_vis","xcorr","mse");
    else
        save("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\participant_id_results_no_filters.mat","centroids","features","tsne_vis","xcorr","mse");
    end
end

