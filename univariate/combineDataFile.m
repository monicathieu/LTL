function t = combineDataFile(file_names, file_path)
t = [];
for i = 1:length(file_names)
    thisFname = file_names{i};
    %load data file
    if isequal(file_names{i}(end-3:end),'.mat')
        load(fullfile(file_path,thisFname));
    else
        error('the files are not in .mat format');
    end
    
    t = [t, theData];
end
       