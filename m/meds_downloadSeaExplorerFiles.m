function meds_downloadSeaExplorerFiles()
if ispc
    roote='w:\';
else
    roote='/u01/rapps/';
end

dire.local=fullfile(roote,'glider_project','data');
dire.remote='glider';

% config.ip='ftp2.dfo-mpo.gc.ca:21';
% config.usr='Glider-ro';
% config.pwd='SVpLzV2YM2';

config.ip='ftp.dfo-mpo.gc.ca';
config.usr='anonymous';
config.pwd='anonymous';

f=ftp(config.ip,config.usr,config.pwd);
cd(f, dire.remote);
cd(f,'realData');
files1=dir(f);
for i=1:length(files1)
    if files1(i).isdir
        cd(f,files1(i).name);
        display(['new directory: ' files1(i).name]);
        files2=dir(f);
        for j=1:length(files2)
            if files2(j).isdir
                cd(f,files2(j).name);
                display(['new directory: ' files2(j).name]);
                files3=dir(f);
                for k=1:length(files3)
                    if files3(k).isdir
                        error('recursion');
                    else
                        f=downfile(f,files3(k).name,fullfile(dire.local,files1(i).name,files2(j).name));
                    end
                end
                display('new directory: ..');
                cd(f,'..');
            else
                f=downfile(f,files2(j).name,fullfile(dire.local,files1(i).name));
               
            end
        end
        cd(f,'..');
        display('new directory: ..');
    else
        f=downfile(f,files1(i).name,fullfile(dire.local));
    end
end
close(f);
exit;
end

function f=downfile(f,fname,localdir)
if ispc
    roote='w:\';
else
    roote='/u01/rapps/';
end

if ~exist(localdir,'dir')
    display (['current at ' localdir]);
    mkdir(localdir);
    system(['chmod 777 ' localdir]);
end
if ~exist(fullfile(localdir,fname), 'file')
    [tr1,tr2,ext]=fileparts(fname);
    if(strcmpi(ext,'.cfg')...
            || strcmpi(ext,'.csv')...
            ||strcmpi(ext,'.log')...
            ||strcmpi(ext,'.msn')...
            || strcmpi(ext,'.kml'))
        ascii(f);
    else
        binary(f);
    end
    mget(f,fname,localdir);
    display(['downloading ' fname]);
    if (strcmpi(ext, '.gz'))
    gunzip ([localdir '/' fname], fullfile(roote,'glider_project','data', 'SeaExplorer', 'input'));
    end
else
    display(['skipping ' fname]);
end
end