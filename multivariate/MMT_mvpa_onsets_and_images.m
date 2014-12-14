function [S] = MMT_mvpa_onsets_and_images(S)

if S.xval
    [onsets_train, S.idxOnsets_train_in_classifier] = MakeOnsets(S, 'train');
    [onsets_test, S.idxOnsets_test_in_classifier] = MakeOnsets(S, 'test');
    
    for k=1:length(onsets_train)
        S.onsets{k} = union(onsets_train{k}, onsets_test{k});
    end
else
    [onsets_train, S.idxOnsets_train_in_classifier] = MakeOnsets(S, 'train');
    [onsets_test, S.idxOnsets_test_in_classifier] = MakeOnsets(S, 'test');
    
    for k=1:length(onsets_train)
        S.onsets{k} = union(onsets_train{k}, S.durTrain + onsets_test{k});
    end
end

S.onsets_train_in_classifier = onsets_train;
S.onsets_test_in_classifier = onsets_test;

end

function [onsets idxOnsets_in_classifier] = MakeOnsets(S, part)

if strcmp(part, 'train')
    if strncmp(S.taskModel, 'study',5) % workaround because the study mvpa onsets are labeled "onsets_mvpa_new" bc we had to strike one rest timepoint
        thisOns = load(fullfile(S.onsetsTrainDir, 'onsets_mvpa_new')); % this loads all the vars into a struct
    else 
        thisOns = load(fullfile(S.onsetsTrainDir, 'onsets_mvpa'));
    end
    thisConds = S.condsTrain;
	thisCondNums = S.condsTrainNums;
    theseFiles = S.filenames_train;
    theseTRWeights = S.TR_weights{1};
    thisIdx = S.idxTr;
elseif strcmp(part, 'test')
        if strncmp(S.taskModel, 'study',5)
        thisOns = load(fullfile(S.onsetsTrainDir, 'onsets_mvpa_new'));
    else 
        thisOns = load(fullfile(S.onsetsTrainDir, 'onsets_mvpa'));
    end
    thisConds = S.condsTest;
	thisCondNums = S.condsTestNums;
    theseFiles = S.filenames_test;
    theseTRWeights = S.TR_weights{2};
    thisIdx = S.idxTe;
end

onsets = cell(size(thisConds));
for i = 1:length(thisConds)
    for j = 1:length(thisConds{i})
        idxThisCond = find(strcmp (thisOns.names, thisConds{i}{j}));
        if ~isempty(idxThisCond)
            time_idx = floor(thisOns.onsets{idxThisCond}/S.TR) + 1;
            enoughTRs_h = (time_idx + length(theseTRWeights)) <= length(theseFiles);
            %enoughTRs{i} = logical(vertcat(enoughTRs{i}, enoughTRs_h));
            theseOnsets = asRow(thisOns.onsets{idxThisCond});
            if strcmp(S.patternType, 'betas')
                onsets{i} = sort(horzcat(onsets{i}, theseOnsets));
            else
                onsets{i} = sort(horzcat(onsets{i}, theseOnsets(enoughTRs_h)));
            end
        end
    end
end

idxOnsets_in_classifier = asRow(ismember(thisIdx.allTrialOns, [onsets{:}]));
%% here, in the magical MMT specific onsets barfy outer, I will specify aggregatedonset conditions. because we didn't analyze responses at scan, all of the behavioral data with which to label the onsets from scratch is sitting in another place
% S.condsTrainNums specifies the columns to include as relevant conditions, since the cell array doesn't have names associated in the struct

% I think onsets is basically the onsets file we have already made, just only including the TRs we care about? formatted as TR and not seconds

% onsets_of_interest = vertcat(horzcat(onsets{thisCondNums{1}}),horzcat(onsets{thisCondNums{2}})); % recreate the onsets except only include from the relevant conditions (concatenating across subconditions if necessary)

end
