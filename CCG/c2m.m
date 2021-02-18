%function yout = c2m(x, els)
function yout = c2m(x, els)

n = length(x);
%y = zeros([sz n]);
for el =els
    sz = size(x{1}{el});
    ndim = length(sz)+1;
    for i=1:n
        addy = shiftdim(x{i}{el},-1);
        if (i==1)
            y = addy;
        else
            y = cat(1, y, addy);
        end
    end
    y = permute(y, [[2:ndim] 1]);
    if length(els)==1
        yout =y;
    else
        yout{end+1} = y;
    end
end