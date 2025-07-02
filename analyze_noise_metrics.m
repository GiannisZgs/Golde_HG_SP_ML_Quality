clear;
close all;
metrics = load("C:\Users\OneDrive\giann\Desktop\ECG HG paper\results_data\metrics_deviation_from_noise.mat");

metrics = metrics.metrics;
top = 5;
struct_fields = fields(metrics);
signals = ["AgCl","HG1","HG2"];
metric_names = ["trend"];%["trend","vlf","low_freq","mid_freq","high_freq","powerline","powerline_harmonic","ultra_high1_freq","ultra_high2_freq"];
for m = 1:length(metric_names) 
    agcl_ch1 = [];    agcl_ch2 = [];    agcl_ch3 = [];
    hg1_ch1 = [];    hg1_ch2 = [];    hg1_ch3 = [];
    hg2_ch1 = [];   hg2_ch2 = [];   hg2_ch3 = [];
    for f = 1:length(struct_fields)
        fieldd = struct_fields(f); fieldd = fieldd{1};
        
        if str2num(fieldd(2:end)) == 11
            fprintf("Skipping participant %s for incomplete data",fieldd)
            agcl_ch1 = [agcl_ch1 -100];
            agcl_ch2 = [agcl_ch2 -100];
            agcl_ch3 = [agcl_ch3 -100];
            hg1_ch1 = [hg1_ch1 -100];
            hg1_ch2 = [hg1_ch2 -100];
            hg1_ch3 = [hg1_ch3 -100];
            hg2_ch1 = [hg2_ch1 -100];
            hg2_ch2 = [hg2_ch2 -100];
            hg2_ch3 = [hg2_ch3 -100];
            continue
        end
    
        for s = 1:length(signals)
            if s==1
                signal_label = 'AgCl';
            elseif s==2
                signal_label = 'HG1';
            elseif s==3
                signal_label = 'HG2';
            end

            for ch = 1:3
                if s == 1 && ch == 1
                    agcl_ch1 = [agcl_ch1 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                elseif s == 1 && ch == 2
                    agcl_ch2 = [agcl_ch2 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                elseif s == 1 && ch == 3
                    agcl_ch3 = [agcl_ch3 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                elseif s == 2 && ch == 1
                    try
                        hg1_ch1 = [hg1_ch1 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                    catch
                        hg1_ch1 = [hg1_ch1 -100];
                    end
                elseif s == 2 && ch == 2
                    try
                        hg1_ch2 = [hg1_ch2 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                    catch
                        hg1_ch2 = [hg1_ch2 -100];
                    end
                elseif s == 2 && ch == 3 
                    try
                        hg1_ch3 = [hg1_ch3 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                    catch
                        hg1_ch3 = [hg1_ch3 -100];
                    end
                elseif s == 3 && ch == 1 
                    hg2_ch1 = [hg2_ch1 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                elseif s == 3 && ch == 2 
                    hg2_ch2 = [hg2_ch2 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                elseif s == 3 && ch == 3 
                    hg2_ch3 = [hg2_ch3 metrics.(fieldd).(signal_label).([char(metric_names(m)),'_deviation']).(['ch',num2str(ch)])];
                end
            end    
        end
    end
end

%% All metrics need to be close to 1
agcl_ch1 = abs(agcl_ch1 - 1);
agcl_ch2 = abs(agcl_ch2 - 1);
agcl_ch3 = abs(agcl_ch3 - 1);
hg1_ch1 = abs(hg1_ch1 - 1);
hg1_ch2 = abs(hg1_ch2 - 1);
hg1_ch3 = abs(hg1_ch3 - 1);
hg2_ch1 = abs(hg2_ch1 - 1);
hg2_ch2 = abs(hg2_ch2 - 1);
hg2_ch3 = abs(hg2_ch3 - 1);

[top_agcl_ch1_values,top_agcl_ch1] = sort(agcl_ch1,'ascend');
top_agcl_ch1 = top_agcl_ch1(1:top); top_agcl_ch1_values = top_agcl_ch1_values(1:top);

[top_agcl_ch2_values,top_agcl_ch2] = sort(agcl_ch2,'ascend');
top_agcl_ch2 = top_agcl_ch2(1:top); top_agcl_ch2_values = top_agcl_ch2_values(1:top);

[top_agcl_ch3_values,top_agcl_ch3] = sort(agcl_ch3,'ascend');
top_agcl_ch3 = top_agcl_ch3(1:top); top_agcl_ch3_values = top_agcl_ch3_values(1:top);

[top_hg1_ch1_values,top_hg1_ch1] = sort(hg1_ch1,'ascend');
top_hg1_ch1 = top_hg1_ch1(1:top); top_hg1_ch1_values = top_hg1_ch1_values(1:top);

[top_hg1_ch2_values,top_hg1_ch2] = sort(hg1_ch2,'ascend');
top_hg1_ch2 = top_hg1_ch2(1:top); top_hg1_ch2_values = top_hg1_ch2_values(1:top);

[top_hg1_ch3_values,top_hg1_ch3] = sort(hg1_ch3,'ascend');
top_hg1_ch3 = top_hg1_ch3(1:top); top_hg1_ch3_values = top_hg1_ch3_values(1:top);

[top_hg2_ch1_values,top_hg2_ch1] = sort(hg2_ch1,'ascend');
top_hg2_ch1 = top_hg2_ch1(1:top); top_hg2_ch1_values = top_hg2_ch1_values(1:top);

[top_hg2_ch2_values,top_hg2_ch2] = sort(hg2_ch2,'ascend');
top_hg2_ch2 = top_hg2_ch2(1:top); top_hg2_ch2_values = top_hg2_ch2_values(1:top);

[top_hg2_ch3_values,top_hg2_ch3] = sort(hg2_ch3,'ascend');
top_hg2_ch3 = top_hg2_ch3(1:top); top_hg2_ch3_values = top_hg2_ch3_values(1:top);

%% By inspecting the each sensor's channel best power ratios, we can see which participants have consistently good ratios
%HG1: 8 (x4), 10(x4), 23(x4), 20 (x2), 1, 6(x2), 39, 40, 15, 25, 3
%HG2: 9(x5), 6(x3), 41(x3), 3(x3), 38,  13(x2),  15,  31, 27, 7, 42
%AgCl: 1(x4), 13(x4), 25(x3), 3, 38(x2), 43, 10,  21, 22, 20, 7, 40, 31, 15

