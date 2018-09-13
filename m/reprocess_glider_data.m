function reprocess_glider_data (a_glider_directory)
if ispc
    roote='w:\';
else
    roote='/u01/rapps/';
end
files = ls (a_glider_directory);
for i = 1: size (files, 1)
    [tr1,tr2,ext]=fileparts(files(i,:));
    if (strcmpi(ext, '.gz'))
    input = [a_glider_directory files(i,:)];
    output = fullfile(roote,'glider_project','data', 'SeaExplorer', 'old');
    fprintf('Processing file : %s\n', input);
    "C:\Program Files\7-Zip\7z.exe" e input -o output
    gunzip (input, fullfile(roote,'glider_project','data', 'SeaExplorer', 'old'));
    end
end

end