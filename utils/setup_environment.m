function setup_environment()
    % Sets up the MATLAB environment for the Golde Hydrogel project
    % Adds all utility directories to the MATLAB path
    % Ensures all required functions and data are accessible
    
    %Get the repository root 
    [repo_root, ~, ~] = fileparts(pwd);
    
    %Add all utility dirs to path
    addpath(genpath(fullfile(repo_root,'utils')));
    addpath(genpath(fullfile(repo_root,'utils','clustering')));
    addpath(genpath(fullfile(repo_root,'utils','ECG_preprocessing')));
    addpath(genpath(fullfile(repo_root,'utils','feature_extraction')));
    addpath(genpath(fullfile(repo_root,'utils','other')));
    addpath(genpath(fullfile(repo_root,'data')));
    
    print_license_info()
    
    fprintf('MATLAB paths set up complete \n')
end

