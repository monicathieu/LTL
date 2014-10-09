function [] = LT_mod_est(par, thisAnalysis)
%LT_mod_est(par, thisAnalysis)
origdir = pwd;
if ~isstruct(par)
    par = LT_Params(par);
end

global defaults
if isempty(defaults)
    spm_defaults;
end
if ~isfield(defaults,'modality')
    defaults.modality = 'FMRI';
end

SPM = [];

thisAnalysisDir = fullfile(par.analysisdir,thisAnalysis);

cd(thisAnalysisDir);
load('SPM.mat');

% get rid of factor field??????
if isfield(SPM,'factor')
SPM = rmfield(SPM,'factor');
end

%% REML Estimation
SPM = spm_spm(SPM);

%% Automatically set up contrasts for factorial designs
if isfield(SPM,'factor')
    if SPM.factor(1).levels > 1
        % don't both if you've only got 1 level and 1 factor
        cons = spm_design_contrasts(SPM);

        %Create F-contrasts
        for i=1:length(cons)
            con  = cons(i).c;
            name = cons(i).name;
            STAT = 'F';
            [c,I,emsg,imsg] = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
            if all(I)
                DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
            else
                DxCon = [];
            end
            if isempty(SPM.xCon),
                SPM.xCon = DxCon;
            else
                SPM.xCon(end+1) = DxCon;
            end
            SPM = spm_contrasts(SPM,length(SPM.xCon));
        end

        %Create t-contrasts
        for i=1:length(cons)
            % Create a t-contrast for each row of each F-contrast
            % The resulting contrast image can be used in a 2nd-level analysis
            Fcon  = cons(i).c;
            nrows = size(Fcon,1);
            STAT  = 'T';
            for r=1:nrows,
                con = Fcon(r,:);
                str = cons(i).name;
                if ~isempty(strmatch('Interaction',str))
                    name = ['Positive ',str,'_',int2str(r)];
                else
                    sp1  = min(find(isspace(str)));
                    name = ['Positive',str(sp1:end),'_',int2str(r)];
                end

                [c,I,emsg,imsg] = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
                if all(I)
                    DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
                else
                    DxCon = [];
                end
                if isempty(SPM.xCon),
                    SPM.xCon = DxCon;
                else
                    SPM.xCon(end+1) = DxCon;
                end
                SPM = spm_contrasts(SPM,length(SPM.xCon));
            end
        end
    end % if SPM.factor(1).levels > 1
end % if isfield(SPM,'factor')

cd(origdir); 
fprintf('Done\n');