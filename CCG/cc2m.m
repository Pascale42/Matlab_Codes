function yout = c2m(x, els)

n = length(x);
%y = zeros([sz n]);
for el =els
    n1=size(x{1}{el},1);
    n2=size(x{1}{el},2);
    for i1=1:n1
        for i2=1:n2
            
            sz = size(x{1}{el});
            ndim = length(sz)+1;
            for i=1:n
                addy = shiftdim(x{i}{el}{i1,i2},-1);
                if (i==1)
                    y{i1,i2} = addy;
                else
                    y{i1,i2} = cat(1, y{i1,i2}, addy);
                end
            end
            y{i1,i2} = permute(y{i1,i2}, [[2:ndim] 1]);
        end
    end
    if length(els)==1
        yout =y;
    else
        yout{end+1} = y;
    end
end