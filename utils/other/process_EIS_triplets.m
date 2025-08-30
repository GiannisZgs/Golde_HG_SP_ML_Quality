function process_EIS_triplets()
    % Processes EIS triplicate measurements
    % Calculates dielectric parameters including permittivity, modulus, and loss tangent
    % Generates visualizations and exports data to Excel/ASCII format
    
    % process_EIS_triplets: Dynamic EIS triplicate analysis for permittivity,
    % modulus, and loss tangent. Saves results to Excel/ASCII and plots subplots.
    clc; clear; close all;
    fprintf('===== EIS Triplicate Processor =====\n');

    % Ask for electrode name
    electrodeName = input('Enter electrode name: ', 's');

    % Select 3 triplicate .txt files
    fprintf('Select 3 triplicate .txt files for %s.\n', electrodeName);
    [files, path] = uigetfile('*.txt', 'Select Triplicate EIS Files', 'MultiSelect', 'on');
    if isequal(files, 0)
        disp('File selection canceled.');
        return;
    end
    if numel(files) ~= 3
        error('Please select exactly 3 files.');
    end

    fullPaths = fullfile(path, files);

    % Import data
    for i = 1:3
        data = readmatrix(fullPaths{i}, 'Delimiter', ';', 'NumHeaderLines', 1);
        freq(:,i) = data(:,2);
        Zr(:,i)   = data(:,3);
        Zi(:,i)   = data(:,4);
    end

    % Verify matching frequencies
    if any(any(abs(diff(freq,1,2))>eps))
        error('Frequency mismatch across triplicates.');
    end
    freq = freq(:,1);

    % Constants
    C0 = 8.854e-12; % F/m
    omega = 2*pi*freq;
    Zmag_sq = Zr.^2 + Zi.^2;

    % Compute dielectric parameters per triplicate
    E_real = Zi ./ (omega .* C0 .* Zmag_sq);
    E_imag = Zr ./ (omega .* C0 .* Zmag_sq);
    M_real = E_real ./ (E_real.^2 + E_imag.^2);
    M_imag = E_imag ./ (E_real.^2 + E_imag.^2);
    tanD   = E_imag ./ E_real;

    % Compute mean & std
    E_real_mean = mean(E_real,2);    E_real_std = std(E_real,0,2);
    E_imag_mean = mean(E_imag,2);    E_imag_std = std(E_imag,0,2);
    M_real_mean = mean(M_real,2);    M_real_std = std(M_real,0,2);
    M_imag_mean = mean(M_imag,2);    M_imag_std = std(M_imag,0,2);
    tanD_mean   = mean(tanD,2);      tanD_std   = std(tanD,0,2);

    % Build output table
    T = table(freq, ...
        E_real_mean, E_real_std, ...
        E_imag_mean, E_imag_std, ...
        M_real_mean, M_real_std, ...
        M_imag_mean, M_imag_std, ...
        tanD_mean, tanD_std, ...
        'VariableNames',{ ...
        'Frequency_Hz', ...
        'E_real_mean','E_real_std', ...
        'E_imag_mean','E_imag_std', ...
        'M_real_mean','M_real_std', ...
        'M_imag_mean','M_imag_std', ...
        'tanD_mean','tanD_std'});

    % Save results
    outExcel = fullfile(path,[electrodeName '_EIS_Params.xlsx']);
    writetable(T,outExcel,'Sheet','Parameters');
    outTxt = fullfile(path,[electrodeName '_EIS_Params.txt']);
    writetable(T,outTxt,'Delimiter','	','FileType','text');
    fprintf('\nSaved: %s\n       %s\n',outExcel,outTxt);

    % Plot parameters in subplots
    params = {'E_real','E_imag','M_real','M_imag','tanD'};
    labels = {'\epsilon''','\epsilon''''','M''','M''''','tan \delta'};
    figure('Name',[electrodeName ' Dielectric'],'Position',[100,100,800,1000]);
    for p=1:length(params)
        subplot(length(params),1,p); hold on;
        y_mean = T.([params{p} '_mean']);
        y_std  = T.([params{p} '_std']);
        errorbar(freq,y_mean,y_std,'o-','LineWidth',1.2);
        set(gca,'XScale','log'); grid on;
        xlabel('Frequency (Hz)'); ylabel(labels{p});
        title([labels{p} ' vs Frequency']);
    end
end
