clear;
close all;

%% Setup environment
[scripts_dir, ~, ~] = fileparts(pwd);
[root_dir, ~, ~] = fileparts(scripts_dir);
setup_script = fullfile(root_dir,'utils','setup_environment.m');
run(setup_script);

%% Implement a clustering analysis to extract average heartbeat profiles
rng(1); % for reproducibility of clustering results
%% Path to the dataset goes here
data_dir = "C:\Users\giann\OneDrive\Desktop\ECG HG paper\results_data\";
save_dir = data_dir;
%% Parameters
use_filters = 1;
plot_cluster_profiles = 0;
plot_tsne = 0;
manually_cleaned = 0;
num_clusters_multiplier = 1; %number of clusters will be the number of participants times this
leads = 3;
sensors = ["AgCl","HG"];
sensors_num = 2;
fs = 200;

if manually_cleaned
    if use_filters
        data = load(fullfile(data_dir,"manually_cleaned_heartbeat_profiles_MA.mat"));
    else
        data = load(fullfile(data_dir,"manually_cleaned_heartbeat_profiles_no_filters.mat"));
    end
else    
    if use_filters
        data = load(fullfile(data_dir,"heartbeat_profiles_MA.mat"));
    else
        data = load(fullfile(data_dir,"heartbeat_profiles_no_filters.mat"));
    end
end

profiles = data.profiling_struct;


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
        [idx_kmeans_pca,centroids_kmeans_pca] = kmeans(reduced_data_99,num_clusters,'Replicates',5,'Start','cluster','MaxIter',1000,'Distance','sqeuclidean');
        reduced_data_with_centroids = [reduced_data_99; centroids_kmeans_pca];

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

        %% Visualize with colored clusters
        if plot_tsne
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
            legend(['Accuracy: ' sprintf('%.2f%%', acc_kmeans * 100)], 'Location', 'northeast');
            suptitle = [char(sensors(s)) ' - Ch.' num2str(ch) ' - Participant Identification - Macro F1: ' sprintf('%.2f%%', macro_f1_kmeans * 100) ' / Accuracy: ' sprintf('%.2f%%', acc_kmeans * 100)];
            sgtitle(suptitle)
        end
        
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
    if plot_cluster_profiles
        figure;
    end
    global_min = 10000000;
    global_max = -10000000;
    for s = 1:2
        for ch = 1:3
            if ~(any(any(isnan(profiles.(participant).(sensors(s)).(['ch',num2str(ch)])))))
                cluster_profile = centroids.(sensors(s)).(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform;
                if plot_cluster_profiles
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

                features.(participant).(sensors(s)).(['ch',num2str(ch)]) = extract_ecg_heartbeat_features(cluster_profile, fs);
            end
            plot_count = plot_count + 1;
            
        end
    end
    if plot_cluster_profiles
        for i = 1:6
            subplot(1,6,i)
            ylim([global_min,global_max])
        end
    end
end


%% Calculate inter-sensor metrics for every channel

for p = 1:length(participants)
    participant = participants(p); participant = participant{1};
    for ch = 1:3
        if ~(any(any(isnan(profiles.(participant).AgCl.(['ch',num2str(ch)])))) || any(any(isnan(profiles.(participant).HG.(['ch',num2str(ch)])))))
            xcorr_cluster = batch_inter_sensor_metrics(centroids.AgCl.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,centroids.HG.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,'xcorr');
     
            mse_cluster = batch_inter_sensor_metrics(centroids.AgCl.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,centroids.HG.(['ch',num2str(ch)]).(participant).kmeans_cluster_mapped_waveform,'nrmse');

            %Save in hierarchical structure
            xcorr.(['ch',num2str(ch)]).(participant) = xcorr_cluster;
            mse.(['ch',num2str(ch)]).(participant) = mse_cluster;
        end
    end
end


%% Save results and data for further processing
if manually_cleaned
    if use_filters
    save(fullfile(save_dir,"manually_cleaned_participant_id_results.mat"),"centroids","features","tsne_vis","xcorr","mse");
else
    save(fullfile(save_dir,"manually_cleaned_participant_id_results_no_filters.mat"),"centroids","features","tsne_vis","xcorr","mse");
    end
else
    if use_filters
        save(fullfile(save_dir,"participant_id_results.mat"),"centroids","features","tsne_vis","xcorr","mse");
    else
        save(fullfile(save_dir,"participant_id_results_no_filters.mat"),"centroids","features","tsne_vis","xcorr","mse");
    end
end

