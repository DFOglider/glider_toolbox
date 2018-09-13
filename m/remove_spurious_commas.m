function sep=remove_spurious_commas(varargin)
%sep=remove_spurious_commas(lin)
%or
%sep=remove_spurious_commas(lin,';')
%or
%sep=remove_spurious_commas(lin,',')
lin=varargin{1};
if nargin==2
    sym=varargin{2};
else
    sym=',';
end
if lin(end)~=sym
    lin(end+1)=sym;
end
sep=[0 find(lin==sym)];
ok=find(lin=='"');
if ~isempty(ok)
    i=find(sep<ok(1),1,'last');
    while i<=length(sep)-1 %read header
        tem=lin(sep(i)+1:sep(i+1)-1);
        if length(tem)>1 && tem(1)=='"' && tem(end)~='"'
            sep=sep([1:i i+2:end]);
        else
            i=i+1;
        end
    end
end
