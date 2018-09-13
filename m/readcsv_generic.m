function s=readcsv_generic(varargin)
%function s=readcsv_generic(filename,nlintoskip,',')
s=[];
infile=varargin{1};
fid=fopen(infile,'r');
if nargin>1
    for i=1:varargin{2} %skip x lines
        fgetl(fid);
    end
end
if nargin>2
    symbo=varargin{3};
else
    symbo=',';
end
lin=[];
while isempty(lin)
    lin=strtrim(fgetl(fid));
end
sep=[0 find(lin==symbo)];
sep=remove_spurious_commas(lin,symbo);
for i=1:length(sep)-1 %read header
    if nargin<3
        tem=lower(strtrim(lin(sep(i)+1:sep(i+1)-1)));
        tem(tem=='$')='S';
        tem(tem=='#')='N';
        tem=makefieldable(tem);
        if ~isempty(tem)
            column{i}=tem;
        else
            column{i}=['c' num2str(i)];
        end
    elseif nargin==3
        column{i}=['c' num2str(i)];
    end    
    column{i}=column{i}(column{i}~='"');
    column{i}=strrep(column{i},'''','');
end
%make sure column names are unique
[ucolumn,i]=unique(column);
j=setdiff(1:length(column),i);
for i=1:length(j)
    ok=strmatch(column(j(i)),column);
    for k=1:length(ok)
        column{ok(k)}=[column{ok(k)} '_' sprintf('%2.2i',k)];
    end
end   
lc=length(column);
j=0;
while ~feof(fid)
    j=j+1;
    if mod(j,1000)==0
        j
    end
    lin=fgetl(fid);
    if ~isempty(lin)
        sep=[0 find(lin==symbo)];
        sep=remove_spurious_commas(lin,symbo);
        for i=1:length(sep)-1
            if i>lc
                column{i}=['c' num2str(i)];
            end
            if diff(sep(i:i+1))==1
                s(j).(column{i})='NaN';
            else
                s(j).(column{i})=lin(sep(i)+1:sep(i+1)-1);
            end
            if length(s(j).(column{i}))>2 && strcmp(s(j).(column{i})([1 end]),'""')
                s(j).(column{i})=s(j).(column{i})(2:end-1);
            end
        end
    end
end
fclose(fid);