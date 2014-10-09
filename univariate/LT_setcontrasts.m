function [] = LT_setcontrasts(par,analysis,overwrite)

origdir = pwd;
if nargin < 3
    overwrite = input('do you want to overwrite existing contrasts? (1 = yes; 0 = no)');
end
if ~isstruct(par) 
    par = LT_Params(par);
end

STAT = par.constat;

thisAnalysisDir = fullfile(par.analysisdir, analysis);
cd(thisAnalysisDir);

fprintf('\nLoading SPM...'); load SPM
fprintf('Done\n');

if overwrite
    SPM.xCon = [];
    delete('con*');
    delete('spmT*');
end


Xsize = size(SPM.xX.xKXs.X,2);

regNames = SPM.xX.name;

idx.bf1 = ~cellfun('isempty', strfind(regNames,'bf(1)'));

switch analysis
    case 'analysisByPerf_noPhase'
        idx.allHit = ~cellfun('isempty', strfind(regNames,'hit'));
        idx.hitSrcCor = ~cellfun('isempty', strfind(regNames,'hitSrcCor'));
        idx.hitSrcInc = ~cellfun('isempty', strfind(regNames,'hitSrcInc'));
        idx.hitNoSrc = ~cellfun('isempty', strfind(regNames,'hitNoSrc'));
        idx.miss = ~cellfun('isempty', strfind(regNames,'miss'));
        idx.CR = ~cellfun('isempty', strfind(regNames,'CR'));
        idx.FA = ~cellfun('isempty', strfind(regNames,'FA'));

        

        idx.respObj = ~cellfun('isempty', strfind(regNames,'respObj'));
        idx.respFace = ~cellfun('isempty', strfind(regNames,'respFace'));
        idx.respScene = ~cellfun('isempty', strfind(regNames,'respScene'));

        idx.obj = ~cellfun('isempty', strfind(regNames,'_obj_'));
        idx.face = ~cellfun('isempty', strfind(regNames,'_face_'));
        idx.scene = ~cellfun('isempty', strfind(regNames,'_scene_'));

        con.all_gr_Baseline = idx.bf1;
        con.hit_gr_CR = idx.bf1 .* balCon_SimpFx(idx.allHit - idx.CR);
        con.hit_gr_Miss = idx.bf1 .* balCon_SimpFx(idx.allHit - idx.miss);
        con.hitSrcHit_gr_hitSrcMiss = idx.bf1 .* balCon_SimpFx(idx.hitSrcCor - idx.hitSrcInc);
        con.hitSrcHit_gr_hitItemHit = idx.bf1 .* balCon_SimpFx(idx.hitSrcCor - idx.hitNoSrc);
    case 'ByPerf_noPhase_simple'
        idx.allHit = ~cellfun('isempty', strfind(regNames,'hit'));
        idx.hitSrcCor = ~cellfun('isempty', strfind(regNames,'hitSrcCor'));
        idx.hitSrcInc = ~cellfun('isempty', strfind(regNames,'hitSrcInc'));
        idx.hitNoSrc = ~cellfun('isempty', strfind(regNames,'hitNoSrc'));
        idx.miss = ~cellfun('isempty', strfind(regNames,'miss'));
        idx.CR = ~cellfun('isempty', strfind(regNames,'CR'));
        idx.FA = ~cellfun('isempty', strfind(regNames,'FA'));

        con.all_gr_Baseline = idx.bf1;
        con.hit_gr_CR = idx.bf1 .* balCon_SimpFx(idx.allHit - idx.CR);
        con.hit_gr_Miss = idx.bf1 .* balCon_SimpFx(idx.allHit - idx.miss);
        con.hitSrcHit_gr_hitSrcMiss = idx.bf1 .* balCon_SimpFx(idx.hitSrcCor - idx.hitSrcInc);
        con.hitSrcHit_gr_hitItemHit = idx.bf1 .* balCon_SimpFx(idx.hitSrcCor - idx.hitNoSrc);
    
    
    case 'analysisByPerfAndPhase'
        idx.hits = ~cellfun('isempty', strfind(regNames,'hit'));
        idx.phase1 = ~cellfun('isempty', strfind(regNames,'Phase1'));
        idx.phase3 = ~cellfun('isempty', strfind(regNames,'Phase3'));
        idx.phase4 = ~cellfun('isempty', strfind(regNames,'Phase4'));
        idx.CR = ~cellfun('isempty', strfind(regNames,'CR'));
        idx.phase3hits = idx.hits .* idx.phase3;
        idx.phase1hits = idx.hits .* idx.phase1;
        idx.phase3and1hits =idx.phase1hits + idx.phase3hits;
        idx.phase4hits = idx.hits .* idx.phase4;
        con.phase4hits_gr_phase3and1hits = idx.bf1 .* balCon_SimpFx(idx.phase4hits - idx.phase3and1hits);
        con.hit_gr_CR = idx.bf1 .* balCon_SimpFx(idx.hits - idx.CR);
end


fn_con = fieldnames(con);

for f = 1:length(fn_con)
    cnames{f} = fn_con{f};
    cvals{f} = con.(fn_con{f});
end

% preallocate
con_name(1:length(cnames)) = {''};
con_vals = cell(1, length(cnames));

for Tt = 1:length(cnames)
    con_name{Tt} = cnames{Tt};

    con_vals{Tt} = double(cvals{Tt});
    
     con_fileName{Tt} = sprintf('con_%04d.img',Tt); 
    
end

%save the contrasts as a structure for future reference
conStruct.con_name = con_name;
conStruct.con_vals = con_vals;
conStruct.con_fileName = con_fileName; 
save conStruct conStruct;

%%put the contrasts into SPM struct & write file
fprintf('\nBeginning contrasts on subject %s\n', par.substr);


cfid = fopen('conlist','wt');
fprintf(cfid, 'Contrasts for Sub %s\nLast run on %s\n', par.substr, date);
fprintf(cfid, 'Regressor Names: %s', strjoin(regNames, ',  '));

% Loop over created contrasts
for k=1:length(con_vals)

    % Basic checking of contrast
    [c,I,emsg,imsg] = spm_conman('ParseCon',con_vals{k},SPM.xX.xKXs,STAT);
    if ~isempty(emsg)
        %disp(emsg);
        warning('Contrast Impossible To Specify. Dummy Contrast Used!');
        [c,I,emsg,imsg] = spm_conman('ParseCon',double(idx.ME),SPM.xX.xKXs,STAT);
        disp(imsg)
    else
        disp(imsg);
    end;

    % Fill-in the contrast structure
    if all(I)
        DxCon = spm_FcUtil('Set',con_name{k},STAT,'c',c,SPM.xX.xKXs);
    else
        DxCon = [];
    end

    % Append to SPM.xCon. SPM will automatically save any contrasts that
    % evaluate successfully.
    if isempty(SPM.xCon)
        SPM.xCon = DxCon;
    elseif ~isempty(DxCon)
        SPM.xCon(end+1) = DxCon;
    end
    SPM = spm_contrasts(SPM,length(SPM.xCon));
    
    lopsided = '';
    if sum(con_vals{k}) ~= 0
        lopsided = 'WARNING: THIS IS A LOPSIDED CONTRAST!';
        fprintf('\n%s\n\n',lopsided);
    end
    fprintf(fopen('conlist','at'),'%d: %s %s\n%s\n\n', k, con_name{k}, lopsided, num2str(con_vals{k}));
end

fclose(cfid);
%copyfile('conlist',[par.logdir filesep 'conlist-' date]);

cd(origdir);
return;

end

function con = balCon_SimpFx(con)
   pos = find(con > 0);
   neg = find(con < 0);
   nPos = length(pos);
   nNeg = length(neg);
   if nPos > nNeg
      con(neg) = con(neg) * nPos/nNeg;
   else
      con(pos) = con(pos) * nNeg/nPos;
   end
end

function con = padConWithZeros( cIn, Xsize )

conLength = length(cIn);
nZeros = Xsize - conLength;
con = [cIn zeros(1,nZeros)];
end