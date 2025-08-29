function metrics = deviation_from_noise_psd(ecg_raw, ecg_processed, fs, target_band, plot_flag,sensor,participant,same_sensor)
% ecg_data: [samples × channels]
% fs: sampling frequency in Hz
% plot_flag: true/false for PSD plots

    if nargin < 4 || isempty(target_band)
        compute_all = true;
    else
        compute_all = false;
    end

    if nargin < 5
        plot_flag = false;
    end

    [L, N] = size(ecg_raw);  % L = samples, N = channels
    
    bands = {
        'trend',        [0.0 0.3];    
        'vlf',          [0.05 0.5];
        'low',          [0.5 5];
        'mid',          [5 15];
        'high',         [15 40];
        'powerline',    [48 52];
        'powerline_harmonic', [78,82];
        'ultra_high1',   [52 78];
        'ultra_high2',   [82 fs/2];
    };
    band_names_plot = ["Baseline Wander","VLF","LF","MF","HF","Powerline","Powerline 1st Harmonic","UHF1","UHF2"];
    band_names = bands(:,1);
    band_ranges = bands(:,2);

    % Select band(s) of interest
    if ~compute_all
        idx = find(strcmp(band_names, target_band));
        if isempty(idx)
            error('Invalid band name. Choose from: %s', strjoin(band_names, ', '));
        end
        band_names = band_names(idx);
        band_ranges = band_ranges(idx);
    end
    
    all_psd_data = struct();
    for ch = 1:N
        x_raw = ecg_raw(:, ch);
        x_proc = ecg_processed(:, ch);
        
        %% Normalize before calculating power
        x_raw = (x_raw - mean(x_raw))/std(x_raw);
        x_proc = (x_proc - mean(x_proc))/std(x_proc);
        
        [Pxx_raw, f] = periodogram(x_raw, [], [], fs);
        [Pxx_proc, ~] = periodogram(x_proc, [], [], fs);

        total_power_raw = trapz(f, Pxx_raw);
        total_power_proc = trapz(f, Pxx_proc);

        metrics.total_power_raw.(['ch',num2str(ch)]) = total_power_raw;
        metrics.total_power_proc.(['ch',num2str(ch)]) = total_power_proc;

        for b = 1:length(band_names)
            band = band_ranges{b};
            idx_f = f >= band(1) & f <= band(2);

            % Compute band power
            power_raw = trapz(f(idx_f), Pxx_raw(idx_f));
            power_proc = trapz(f(idx_f), Pxx_proc(idx_f));

            % Store raw, processed, and ratio
            metrics.([band_names{b} '_raw']).(['ch',num2str(ch)]) = power_raw;
            metrics.([band_names{b} '_proc']).(['ch',num2str(ch)]) = power_proc;
            metrics.([band_names{b} '_reduction_ratio']).(['ch',num2str(ch)]) = power_proc / (power_raw + 1e-12);
        end

        % Optional plot
        if plot_flag
            
        
            fig = figure('Name', ['PSD - Participant ' participant ' - Sensor ' sensor ', Channel ' num2str(ch)], 'Position', [100 100 1200 600]);
            sgtitle(['Full PSD - Participant ' participant ' - Sensor ' sensor ', Channel ' num2str(ch)]);
            % Full PSD plot
            subplot(2, 5, 1);  % left-most plot
            if same_sensor
                plot(f, 10*log10(Pxx_raw), 'r', 'DisplayName', 'Raw'); hold on;
                plot(f, 10*log10(Pxx_proc), 'b', 'DisplayName', 'Processed');
            else
                plot(f, 10*log10(Pxx_raw), 'r', 'DisplayName', 'AgCl Raw'); hold on;
                plot(f, 10*log10(Pxx_proc), 'b', 'DisplayName', 'HG Raw');
            end
            xlabel('Frequency (Hz)');
            ylabel('PSD (dB/Hz)');
            ylim([-130,10])
            title('Full PSD');
            grid on;
            legend('show','Location','southeast');

            % Annotate bands with xlines and center text
%             y_lims = ylim;
%             for b = 1:length(band_names)
%                 band = band_ranges{b};
%                 xline(band(1), '--k', 'HandleVisibility', 'off');
%                 xline(band(2), '--k', 'HandleVisibility', 'off');
%                 center_freq = mean(band);
%                 text(center_freq, y_lims(2)*0.5, band_names{b}, ...
%                     'HorizontalAlignment', 'center', ...
%                     'FontSize', 9, ...
%                     'BackgroundColor', 'white', ...
%                     'EdgeColor', 'black');
%             end

            % Threshold for zoom-in band width
            %narrow_band_thresh = 5; % Hz

            % Add zoomed-in subplots for narrow bands
            zoom_plot_idx = 2;  % start from subplot 2
            for b = 1:length(band_names)
                band = band_ranges{b};
                band_width = band(2) - band(1);

                subplot(2, 5, zoom_plot_idx);
                idx = f >= band(1) & f <= band(2);
                plot(f(idx), 10*log10(Pxx_raw(idx)), 'r'); hold on;
                plot(f(idx), 10*log10(Pxx_proc(idx)), 'b');
                xlabel('Frequency (Hz)');
                ylabel('PSD (dB/Hz)');
                ylim([-130,10])
                xlim([band(1) band(2)])
                title([char(band_names_plot(b)) ' band (' num2str(band(1)) '-' num2str(band(2)) ' Hz)']);
                grid on;
                zoom_plot_idx = zoom_plot_idx + 1;
            end
            if same_sensor
                filename = sprintf('C:\\Users\\giann\\OneDrive\\Desktop\\ECG HG paper\\matlab_plots\\manually_cleaned_power_spectral_density_same_sensor\\PSD_%s_%s_ch_%d.png', participant, sensor, ch);
            else
                filename = sprintf('C:\\Users\\giann\\OneDrive\\Desktop\\ECG HG paper\\matlab_plots\\manually_cleaned_power_spectral_density_diff_sensor\\PSD_%s_%s_ch_%d.png', participant, sensor, ch);
            end
            saveas(fig, filename);  
            
            %% Save data used for this plot
            channel_key = sprintf('channel_%d', ch);
            psd_data = struct();
            psd_data.channel = ch;
            psd_data.freq = f;
            psd_data.Pxx_raw = Pxx_raw;
            psd_data.Pxx_proc = Pxx_proc;
            psd_data.fs = fs;

            % Store band powers
            for b = 1:length(band_names)
                band = band_ranges{b};
                idx_f = f >= band(1) & f <= band(2);
                psd_data.bands.(band_names{b}).range = band;
                psd_data.bands.(band_names{b}).power_raw = trapz(f(idx_f), Pxx_raw(idx_f));
                psd_data.bands.(band_names{b}).power_proc = trapz(f(idx_f), Pxx_proc(idx_f));
            end

            all_psd_data.(channel_key) = psd_data;

        end
    end  
    if plot_flag
        json_str = jsonencode(all_psd_data);
        if same_sensor
            filename = sprintf('C:\\Users\\giann\\OneDrive\\Desktop\\ECG HG paper\\matlab_plots\\manually_cleaned_power_spectral_density_same_sensor_data\\PSD_%s_%s.json', participant, sensor);
        else
            filename = sprintf('C:\\Users\\giann\\OneDrive\\Desktop\\ECG HG paper\\matlab_plots\\manually_cleaned_power_spectral_density_diff_sensor_data\\PSD_%s_%s.json', participant, sensor);
        end
        fid = fopen(filename, 'w');
        fwrite(fid, json_str, 'char');
        fclose(fid);
    end
end