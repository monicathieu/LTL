function [S idxTr idxTe par]= MMT_mvpa_params(subj_id, task)

% establish parameters for mvpa analysis
% <subj_id> - identifier for the given subject. Can be numerical or a
% string


%% establish general parameters
idxTr = [];
idxTe = [];
par = LT_Params(subj_id, task);
S.exp_name = 'MMT';

%% directories
S.subj_id = par.substr;
S.expt_dir = '/biac4/wagner/wagner/mthieu/MMT_fMRI/fmri_data/';
S.mvpa_dir = [S.expt_dir  S.subj_id '/mvpa_results'];
S.anat_dir = [S.expt_dir  S.subj_id '/anatomicals_study']; % if using within task?
S.importance_maps_dir=[S.expt_dir 'mvpa_results/ImpMaps_' date  ];
S.group_mvpa_dir = [S.expt_dir 'mvpa_results'];


S.workspace_dir = [par.subdir '/' 'mvpa_workspace']; % check me

%% preprocessing
S.preprocType = 'spm'; % 'spm' for spm preprocessing, 'knk' for kendrick preprocessing
S.loadAndPreprocFunctName = 'JR_mvpa_load_and_preprocess_raw_data'; % don't need to use this hopefully
%% tasks
S.trainTask = 'study'; % so here is where you would specify train on localizer
S.testTask = 'study';

%% cross-validation scheme
if strcmp(S.trainTask, S.testTask)
    S.xval = 1;
    S.thisSelector =  'randomNFold_xval'; % cross validation. 
else
    S.xval = 0; % AHA great this is where the across task is specified???
    S.thisSelector = 'TrainTestOneIterGroup'; % train on one group, test on another
end

%% information specific to the training and testing of specific sets of data

% parTr = par for the training data
% S.onsetsTrainDir = location of training onsets
% S.condsTrain = conditions on which to train
% S.dnCondsTrain = conditions which which to denoise, if denoising is used
% S.TrainRuns = runs of data on which to train
% S.durTrain = duration of training
% S.filenames_train = names of data images to use for training
% idxTr = behavioral information for training task
    
%% training
% every train "task" is actually a specific condition combination of interest WITHIN said task
switch S.trainTask
    case 'loc_face_vs_obj'
        parTr = LT_Params(subj_id, 'loc'); % what is parTr?
	S.taskModel = 'loc';
        S.onsetsTrainDir = [S.expt_dir S.subj_id '/results_localizer/'];
	S.condsTrain = {{'Faces'}  {'Objs'}} ;
	S.condsTrainNums = {3 2}; % these must be relative to the onsets file they come from 
        S.TrainRuns = par.scansSelect.(par.task);
        S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
        S.filenames_train = vertcat(par.swascanfilesByRun.lab{S.TrainRuns});
        [~, ~,idxTr] = LTLab_fMRIBehAnalysis_Ret(par);
    case 'study_face_vs_obj_all'
        parTr = LT_Params(subj_id, 'study');
        S.onsetsTrainDir = [S.expt_dir S.subj_id '/results01/'];
	S.taskModel = 'study1';
	S.condsTrain = {{'Faces_All'} {'Objs_All'}};
	S.condsTrainNums = {[3,4] [5,6]};
        S.TrainRuns = par.scansSelect.(par.task);
        S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
        S.filenames_train = vertcat(par.swascanfilesByRun.lab{S.TrainRuns});
        [~, ~,idxTr] = LTLab_fMRIBehAnalysis_Ret(par);

%% testing
switch S.testTask
    case 'study_face_vs_obj_all'
	parTe = LT_Params(subj_id, 'study');
	S.taskModel = 'study1';
        S.onsetsTestDir = [S.expt_dir S.subj_id '/results01/'];
        S.condsTest = {{'Face_All'}  {'Obj_All'} } ;
	S.condsTestNums = {[3,4] [5,6]};
        S.TestRuns = par.scansSelect.(par.task);
        S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
        S.filenames_test = vertcat(par.swascanfilesByRun.lab{S.TestRuns});
        % %     S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
        [~, ~,idxTe] = LTLab_fMRIBehAnalysis_Ret(par);
    case 'study_face_vs_obj_hits'
        parTe = LT_Params(subj_id, 'study');
	S.taskModel = 'study1';
        S.onsetsTestDir = [S.expt_dir S.subj_id '/results01/'];
        S.condsTest = {{'Face_Hits'}  {'Obj_Hits'}} ;
	S.condsTestNums = {3 5};
        S.TestRuns = par.scansSelect.(par.task);
        S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
        S.filenames_test = vertcat(par.swascanfilesByRun.lab{S.TestRuns});
        [~, ~,idxTe] = LTLab_fMRIBehAnalysis_Ret(par);
    case 'study_face_vs_obj_misses'
        parTe = LT_Params(subj_id, 'study');
	S.taskModel = 'study1';
        S.onsetsTestDir = [S.expt_dir S.subj_id '/results01/'];
        S.condsTest = {{'Face_Misses'}  {'Obj_Misses'}} ;
	S.condsTestNums = {4 6};
        S.TestRuns = par.scansSelect.(par.task);
        S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
        S.filenames_test = vertcat(par.swascanfilesByRun.lab{S.TestRuns});
        [~, ~,idxTe] = LTLab_fMRIBehAnalysis_Ret(par);
end

S.condnames = S.condsTrain;
S.regName = 'conds';


%% Smoothing Parameters
S.funcType = 2; % leave this alooonnnneeeeee
S.smoothTxt = { 'unsmoothed' 'smoothed' 'native'};
switch S.funcType
    case 1
    par.filesForPatterns = par.wascanfiles;
    case 2
    par.filesForPatterns = par.swascanfiles;
    case 3
    par.filesForPatterns = par.ascanfiles;     
end

%% specify which files to load for classification
if S.xval
    S.filenames = S.filenames_train;
else
%     S.filenames_h{1} = S.filenames_train;
%     S.filenames_h{2} = S.filenames_test;
    S.filenames = [S.filenames_train; S.filenames_test];
end
    S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);

%% Runs Parameters

%S.runs_vector - number of volumes per each run
if S.xval
    S.runs_vector =  [par.(par.task).numvols(S.TrainRuns)];
else
    S.runs_vector =  [parTr.(parTr.task).numvols(S.TrainRuns) par.(par.task).numvols(S.TestRuns)];
    S.TrainTestOrder = [ones(size(S.TrainRuns)) 2*ones(size(S.TestRuns))];
end


S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = par.TR;

%% Volume Parameters
S.vol_info = spm_vol(fullfile(par.funcdir, 'study01', 'swavol0005.nii')); %get functional data resolution info for spm .img writing

% % S.roiWithNonTaksVoxels = fullfile(par.anatdir, 'tnativeOccTemp.nii');
% % S.roiWithNonTaksVoxelsName = 'tnativeOccTemp.nii';
% can I change this so I can specify the mask type earlier?
S.roi_file = [S.expt_dir '/Masks/rSEPT09_MVPA_MASK_resliced4mm.nii'];
S.roi_name = 'rSEPT09_MVPA_MASK_resliced4mm.nii';

S.secondaryMask = []; % secondary mask (the specific classification mask)
%masks the primary data loaded in the workspace. [] = no secondary mask.
%S.secondaryMask = [S.expt_dir 'rLTL_loc_Faces_gr_Scenes05k5.img'];

%% Workspace Parameters
S.use_premade_workspace = 1;
% S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.funcType} '_train_' S.trainTask '_test_' S.testTask S.preprocType '.mat']);
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.funcType} '.mat']);



%% Pattern names
S.patternType = 'raw'; %'raw' or 'betas'
% S.preprocPatName = 'spiral_dn';
% S.preprocPatName = 'patsAllVox_z_dn';
S.preprocPatName = 'epi_d_hp_z';
S.preprocPatCondensedName = [S.preprocPatName '_condensed'];

if isempty(S.secondaryMask)
    S.preprocPatNameFinalMask = S.preprocPatName;
else
    S.preprocPatNameFinalMask = [S.preprocPatName '_masked'];
end

%% Artifacts
% % S.artFile = (fullfile(par.artrepdir, ['art_global_modified_' par.substr])); %directory where artifact information is contained
S.inactivateArtifacts = 0; %remove artifact trials?

%% Iteration Parameters
S.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
S.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data

%% Balancing Parameters
S.equate_number_of_trials_in_groups = 1; % equate number of trials in conditions 1 and 2
S.numBalancedIts = 1; % number of iterations to run, with different randomization for the balancing

%% Z-Scoring and outlier detection
S.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification 
S.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
S.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
S.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
S.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)

%% Importance Maps
S.generate_importance_maps = 1; %visualize classifier weights
S.generateBetaMaps = 1; %use betas, instead of importance values
S.impType = {'pos' 'neg' 'both' 'raw'}; %importance map types
S.regNames = S.condsTrain;

%% Special types of analysis
S.searchlightAnalysis = 0; % run a searchlight analysis
S.linReg = 0; % run an analysis with a continuous outcome variable

%% Subsample
S.subsampleToMatch = 0; %subsample trials to match quantities across them. is this not a duplication of the balancing flag above?
S.balanceHiAndLowConf = 0;% match N of hi and low confidence trials?

%% voxel interactions
S.includeVoxelInteractions = 0; %include interactions among voxels?   
S.interactionType = 2;
S.intEst = 1;
S.intConcat = 0;
S.intPThresh = .001;
S.intReportIncrement = 100;
S.intFilePrefix = 'intVox';
S.intMaskName = 'interactionMask';
S.intPatName = 'interactions';
S.intGroupName = 'interactionsGroup';
S.intUseIntsWithUniqueInfo = 1;

%% Signal intensity analysis
S.thisSigIntenseSelector = 'randomNFold_xval'; %which selector to use for signal intensity analysis
S.zscoreIntensityVals = 1; % zscore the intensity values?

%% Denoising
S.denoise = 0; %undergo denoising?
S.denoiseOpt.denoisespec = '10001'; %which parts of the glm output do we want to save?

%% Mean Signal Extraction Params
% parameters for selecting the mean signal from a class-specific ROI for each pattern.

S.extractMeanSignal = 0; %1 - do signal extraction. 0 - don't do this. 
S.defineROIsFromANOVAFS = 0; % define ROIs using ANOVA-based feature selection, instead of pre-defining them. 
S.logreg_2Features = 0; %perform a logistic regression, using the two extracted intensity vectors

% % S.ROI1PatName = [S.preprocPatCondensedName '_ROI1'];
% % S.ROI1_name = ['occipitoTemporal_faceVsScene_500vox.img'];
% % S.ROI1_file  = [par.subdir '/analysis_loc_mnem/' S.ROI1_name];
% % 
% % S.ROI2PatName = [S.preprocPatCondensedName '_ROI2'];
% % S.ROI2_name  = ['occipitoTemporal_sceneVsFace_500vox.img'];
% % S.ROI2_file   = [par.subdir '/analysis_loc_mnem/' S.ROI2_name];

%% TR Weighting
%which post-stimulus TRs should be used (and if more than one, averaged
%across) before feeding data to the classifier?  S.TR_weights_set{1} -
%training weights, S.TR_weights_set{2} = testing weights
S.TR_weights_set = {[0 0 .33 .34 .33] [0 0 .33 .34 .33]}; % both train and test are set to late TR (average over the last 3)

%% classifier parameters
S.class_args.train_funct_name = 'train_pLR'; %training function
S.class_args.test_funct_name = 'test_pLR'; %testing function
S.class_args.classType = 'pLR';
S.perfmet_functs = 'perfmet_maxclass'; % performance metric
S.statmap_funct = 'AG_statmap_anova'; % you don't have this function! 
S.nPlsCompsSet = 0; % number of pls components to include. 0 = do not use pls components.
S.nFolds = 5; % number of cross validation iterations

S.class_args.nVox = 0; % number of voxels to select with feature selection e.g. [1000 5000 10000]
S.class_args.libLin = '-q -s 0 -B 1'; %arguments for liblinear
S.class_args.libsvm = '-q -s 0 -t 2 -d 3'; % arguments for libsvm
S.class_args.constant = true; % include a constant term?
S.class_args.prefitWeights = true; 
S.class_args.chooseOptimalPenalty = 0; % cycle through cost parameters in the training set, and chose the optimal one?
S.class_args.penaltyRange = [.001 .005 .01 .05 .1 .5 1 5 10 50 100 500 1000 50000]; % cost parameters to cycle through
S.class_args.radialBasisSelection = [];%[.00001 .0001 .001 .01 .1 1 10];
S.class_args.nFoldsPenaltySelection = 5; % number of cross validation folds for penalty parameter selection. 

S.class_args.penalty = 10; % default penalty if chooseOptimalPenalty=0

end
