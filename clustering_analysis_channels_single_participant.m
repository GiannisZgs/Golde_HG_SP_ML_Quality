clear;
close all;
%% Implement a clustering analysis to extract average heartbeat profiles
%% This should be done across patients, across sensors, across channels
use_filters = 1;
if use_filters
    data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\single_participant_heartbeat_profiles_MA.mat");
else
    data = load("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\single_participant_heartbeat_profiles_no_filters.mat");
end

fs = 200;
profiles = data.profiling_struct;
num_clusters_multiplier = 1; %number of clusters will be the number of participants times this
leads = 3;
sensors = ["AgCl","HG"];
sensors_num = 2;
rng(1); % for reproducibility

%% Take all channels, sensors of all participants, cluster them into N_sensors*N_channels clusters
%Gather all data, assign labels according to channel-sensor

channel_mapping = ["Lead1","Lead2","Lead3"];
for s = 1:sensors_num
    channel = [];
    labels = [];
    for ch=1:leads
        participants = fields(profiles);
        participants_used = 0;
        for p = 1:2 %length(participants)
            participant = participants(p); participant = participant{1};
            if ~any(any(isnan(profiles.(participant).(sensors(s)).(['ch',num2str(ch)]))))
                signal = profiles.(participant).(sensors(s)).(['ch',num2str(ch)]);
                channel = [channel signal];
                labels = [labels; ch*ones(size(signal,2),1)];
            end
        end   
        %centroids.(sensors(s)).(['ch',num2str(ch)]).direct_centroid_waveform = mean(channel,2);
    end

    channel = channel';
    %% Reduce data to visualize - PCA
    %Ensure kept components explain at least 95% of the variability
    %Don't go over 50 components
    [coeff, score, latent, ~, explained,mu] = pca(channel);
    cumulative_explained = cumsum(explained);
    num_components = find(cumulative_explained >= 99, 1, 'first');
    num_components = min(num_components, 50);
    reduced_data_99 = score(:,1:num_components);

    num_clusters = leads*num_clusters_multiplier;
    %[idx_kmeans_pca,centroids_kmeans_pca] = kmeans(reduced_data_99,num_clusters,'MaxIter',1000); %'Replicates',5,'Start','cluster','MaxIter',1000,'Distance','sqeuclidean');
    [idx_kmeans_pca,centroids_kmeans_pca] = kmeans(reduced_data_99,num_clusters,'Replicates',5,'Start','cluster','MaxIter',1000,'Distance','sqeuclidean');
    reduced_data_with_centroids = [reduced_data_99; centroids_kmeans_pca];

    %% Reduce further to 2 dimensions for visualizing
    combined_tsne = tsne(reduced_data_with_centroids, 'NumDimensions', 2, 'Perplexity', 30);

    tsne_data = combined_tsne(1:end-num_clusters, :);
    tsne_centroids = combined_tsne(end-num_clusters+1:end, :);

    %% Inverse transform of centroids from principal component space to time domain space
    centroids_waveforms = centroids_kmeans_pca * coeff(:,1:num_components)' + mu;

    %% Compute clustering accuracy and map clusters to true labels - Assignments(j) is the true label corresponding to cluster j
    [acc_kmeans,macro_f1_kmeans,micro_f1_kmeans,cluster_assignments] = greedy_cluster_mapping(labels, idx_kmeans_pca);
    for j = 1:length(cluster_assignments)
        centroids.(sensors(s)).(['ch',num2str(j)]).kmeans_cluster_mapped_waveform = centroids_waveforms(j,:)';
    end
    fprintf('%s Unsupervised Channel Identification Accuracy (K-Means): %.2f%%\n', sensors(s), acc_kmeans * 100);
    fprintf('%s Unsupervised Channel Identification Macro-F1 (K-Means): %.2f%%\n', sensors(s), macro_f1_kmeans * 100);
    fprintf('%s Unsupervised Channel Identification Micro-F1 (K-Means): %.2f%%\n', sensors(s), micro_f1_kmeans * 100);

%     channel_identification.(sensors(s)).acc = acc_kmeans;
%     channel_identification.(sensors(s)).macro_f1 = macro_f1_kmeans;
%     channel_identification.(sensors(s)).micro_f1 = micro_f1_kmeans;
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
    %lgd = findobj('type','legend');
    %delete(lgd)
    %legend(['Accuracy: ' sprintf('%.2f%%', acc_kmeans * 100)], 'Location', 'northeast');
    suptitle = [char(sensors(s)) ' - Channel Identification - Macro F1: ' sprintf('%.2f%%', macro_f1_kmeans * 100) ' / Accuracy: ' sprintf('%.2f%%', acc_kmeans * 100)];
    sgtitle(suptitle)
    
    %% Store TSNE data used in plots for further processing
    tsne_vis.(sensors(s)).projection = tsne_data;
    tsne_vis.(sensors(s)).centroids = tsne_centroids;
    tsne_vis.(sensors(s)).kmeans_clusters = idx_kmeans_pca;
    tsne_vis.(sensors(s)).ground_truth = labels;
    tsne_vis.(sensors(s)).accuracy = acc_kmeans;
    tsne_vis.(sensors(s)).macro_f1 = macro_f1_kmeans;
    tsne_vis.(sensors(s)).micro_f1 = micro_f1_kmeans;
end

%% Visualize channel profiles
%Cluster profiles
plot_count = 1;
figure;
global_min = 10000000;
global_max = -10000000;
for s = 1:2
    for ch = 1:3
        cluster_profile = centroids.(sensors(s)).(['ch',num2str(ch)]).kmeans_cluster_mapped_waveform;
        %subplot(length(plot_participants),6,plot_count)
        subplot(1,6,plot_count)
        plot(cluster_profile,'k','LineWidth',1.5)
        xlabel((['ch',num2str(ch)])); %
        if ch == 2
            title(sensors(s))
        end
        if min(cluster_profile) < global_min
            global_min = min(cluster_profile);
        end
        if max(cluster_profile) > global_max
            global_max = max(cluster_profile);
        end
        %subplot(length(plot_participants),6,plot_count)
        %subplot(1,6,plot_count)
        plot_count = plot_count + 1;
        features.(sensors(s)).(['ch',num2str(ch)]) = extract_ecg_heartbeat_features(cluster_profile, fs);
    end
end
for i = 1:6
    subplot(1,6,i)
    ylim([global_min,global_max])
end

%% Calculate inter-sensor metrics for every channel

for ch = 1:3
    cluster_profile_1 = centroids.AgCl.(['ch',num2str(ch)]).kmeans_cluster_mapped_waveform;
    cluster_profile_2 = centroids.HG.(['ch',num2str(ch)]).kmeans_cluster_mapped_waveform;

    xcorr_cluster = batch_inter_sensor_metrics(cluster_profile_1,cluster_profile_2,'xcorr');
    mse_cluster = batch_inter_sensor_metrics(cluster_profile_1,cluster_profile_2,'nrmse');
        
    %Save in hierarchical structure
    xcorr.(['ch',num2str(ch)]) = xcorr_cluster;
    mse.(['ch',num2str(ch)]) = mse_cluster;
end

save("C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\single_participant_channel_id_results.mat","centroids","features","tsne_vis","xcorr","mse");

%sgtitle('Centroid Profiles obtained from Clustering Analysis for 4 participants')

