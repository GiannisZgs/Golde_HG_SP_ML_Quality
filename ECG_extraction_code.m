% you need to add the functions to the relative path also explain the xlsx
% function and the writematrix funtion
% call them if needed

%add ask me toenter the file name instead of changing code everytime
%add more order in where i.e folder
clc;clear
%to name each array [A,B,C]= function('Filename.format)
%to hide its output in script windows you added a semicolon
[header, signal,annotations]=ConvertKHA2MAT('C:\Users\nazmi\OneDrive\Desktop\KMTPatient_00.kha');
% to save that as mat file
save("ECG.mat","header","signal","annotations");
%%load('AimbuECG.mat')
 %%plot(signal(1).data(:)) signal is a structured array of the mat file
 %extracted from the kha file format.
 % signal(row number).data(rangemin:rangemax);
 %signal(1).transducer='TRG' name the column as TRG
 %to verify the correctness of the measurement you can 
 % use smp_in_file ... to check of the no of entries calc with sampling
 % freq
 %header and annot also the patient name and DOB must match the outer txt
 %files

%Displays the all perfect

%ECG LEAD 1
% Load the ECG signal from the .mat file
load('ECG.mat', 'signal');

% Extract the ECG lead data 
ecg_lead_data = signal(1).data;

% Read Amplitude 
amplitude_uV = ecg_lead_data(:);  % Convert to a single column
% Create a time array based on the sampling rate (200 Hz)
sampling_rate = 200;
time = (0:1:length(amplitude_uV)-1) / sampling_rate; % Adjusted length

% Ensure both arrays have the same number of rows
min_length = min(length(time), length(amplitude_uV));
time = time(1:min_length);
amplitude_uV = amplitude_uV(1:min_length);
timep=time';
% Create a table with only the amplitude data
T = table(timep,amplitude_uV, 'VariableNames', {'Time(s)','ECG_LEAD_I'});
% Write the table to an Excel file with the variable name as the header
writetable(T, 'Patient0_agcl.xlsx', 'sheet','ECG lead 1');
disp(['ECG_LEAD I data written to xlsx successfully.']);

%ECG LEAD 2
% Extract the ECG lead data (assuming it's stored in the first column)
ecg_lead_data = signal(2).data;

% Convert amplitude to microvolts (uV)
amplitude_uV = ecg_lead_data(:);  % Convert to a single column

% Ensure both arrays have the same number of rows
min_length = min(length(time), length(amplitude_uV));
time = time(1:min_length);
amplitude_uV = amplitude_uV(1:min_length);
timep=time';
% Create a table with only the amplitude data
T = table(timep,amplitude_uV, 'VariableNames', {'Time(s)','ECG_LEAD_II'});
% Write the table to an Excel file with the variable name as the header
writetable(T, 'Patient0_agcl.xlsx', 'sheet','ECG lead 2');
disp('ECG_LEAD II data written to xlsx successfully.');

%ECG LEAD III
% Extract the ECG lead data (assuming it's stored in the first column)
ecg_lead_data = signal(3).data;

% Convert amplitude to microvolts (uV)
amplitude_uV = ecg_lead_data(:);  % Convert to a single column

% Ensure both arrays have the same number of rows
min_length = min(length(time), length(amplitude_uV));
time = time(1:min_length);
amplitude_uV = amplitude_uV(1:min_length);
timep=time';
% Create a table with only the amplitude data
T = table(timep,amplitude_uV, 'VariableNames', {'time(s)','ECG_LEAD_III'});
% Write the table to an Excel file with the variable name as the header
writetable(T, 'Patient0_agcl.xlsx', 'sheet','ECG lead 3');
disp('ECG LEAD III data written to xlsx successfully.');

% % %  data = load('AimbuECG.mat','signal');
% % % f = fieldnames(data);
% % % time = (0:1:length(data.(f{1}))-1) / 200; % adjust for 200Hz sampling rate
% % % amplitude = data.(f{1}) * 1e6; % convert to uV
% % % 
% % % % Create a cell array with 'ECG LEAD I' as the first row and the data as the following rows
% % % C = [{'ECG LEAD I'}; num2cell([time', amplitude'])];
% % % 
% % % % Write the cell array to an Excel file
% % % writecell(C, 'AimbuECG.xlsx');

%%%%%%%%%%%%%%%%%%
% % Load the ECG signal from the .mat file
% load('AimbuECG.mat', 'signal');
% 
% % Extract the ECG lead (assuming it's stored in the first column)
% ecg_lead = signal.data;
% 
% % Create a time array based on the sampling rate (200 Hz)
% sampling_rate = 200;
% time = (0:1:length(ecg_lead)-1) / sampling_rate;
% 
% % Convert amplitude to microvolts (uV)
% amplitude_uV = ecg_lead;
% 
% % Create a table with time and amplitude
% T = table(time', amplitude_uV, 'VariableNames', {'Time_s', 'ECG_LEAD_I'});
% 
% % Write the table to an Excel file with the variable names as headers
% writetable(T, 'AimbuECG.xlsx');
% % %This Code Get you the ECG data in one column bas only amplitude
% % 
% % % Load the ECG signal from the .mat file
% % load('AimbuECG.mat', 'signal');
% % 
% % % Create a time array based on the sampling rate (200 Hz)
% % sampling_rate = 200;
% % time = (0:1:length(ecg_lead)-1) / sampling_rate;
% % 
% % % Extract the ECG lead data (assuming it's stored in the first column)
% % ecg_lead_data = signal.data;
% % 
% % % Convert amplitude to microvolts (uV)
% % amplitude_uV = ecg_lead_data(:);  % Convert to a single column
% % 
% % % Create a table with only the amplitude data
% % T = table(time,amplitude_uV, 'VariableNames', {'Time_s','ECG_LEAD_I'});
% % 
% % % Write the table to an Excel file with the variable name as the header
% % writetable(T, 'AimbuECG.xlsx', 'Sheet', 1);
% % 
% % disp('Amplitude data written to AimbuECG.xlsx successfully.');

