function [subj] = LTC_mvpa_load_and_preprocess_raw_data_JR(S)
    
    trial_idx = 1;
    for r = 1:length(S.runs_vector)
        runs(trial_idx:trial_idx+S.runs_vector(r)-1)=r;
        trial_idx = trial_idx+S.runs_vector(r);
    end
    S.runs = runs;

    % initialize subj structure
    subj = init_subj(S.exp_name,S.subj_id);

    % load mask file
    subj = load_spm_mask(subj,S.roi_name,S.roi_file);
    
    % load functional data
    subj = load_analyze_pattern(subj,'epi',S.roi_name,S.filenames,'single',true); %use single precision format to save RAM

    % move pattern to hard drive to save RAM (optional)
    %subj = move_pattern_to_hd(subj, 'epi');

    % make runs vector
    subj = init_object(subj,'selector','runs');
    subj = set_mat(subj,'selector','runs',runs);

    % detrend the timeseries data
    subj = detrend_runs(subj,'epi','runs');  % not in mvpa tutorial, but seems important to do

    % move pattern to hard drive to save RAM (optional)
    %subj = move_pattern_to_hd(subj, 'epi_d');
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','epi');
    
    % high-pass filter the timeseries data
    subj = hpfilter_runs(subj,'epi_d','runs',128,2); % remove frequencies below .01 Hz (adjust as desired)

     % clean up workspace
    subj = remove_mat(subj,'pattern','epi_d');
    
    % move pattern to hard drive to save RAM (optional)
    %subj = move_pattern_to_hd(subj, 'epi_d_hp');
    
    % zscore the data from each run
    subj = zscore_runs(subj,'epi_d_hp','runs'); % gives each voxel a mean of 1 and variance of 0 across all timepoints of each run
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','epi_d_hp');
    
    %save final pattern in single precision form (8 sig figs) to save RAM and HD space    
    subj.patterns{end}.mat = single(subj.patterns{end}.mat);
    
    %save the workspace
    if ~exist(S.workspace_dir)
        mkdir(S.workspace_dir);
    end
    cd (S.workspace_dir);
    save (S.workspace);