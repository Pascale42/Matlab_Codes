function fieldsz =  FieldSizes(ins,f);
% if nargin<2 | isempty(f)
%     f = fieldnames(ins(1));
fieldsz=[];
for k=1:length(ins(:))
    v=getfield(ins, {k}, f);
%    if ~isempty(v)
        fieldsz(k,:) = size(v);
%    end
end