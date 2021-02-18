function xi = MakeUniformDistrInv(x,xu,xiu)
xi = zeros(size(xiu));
% [xu, si]= sort(xu);
% x = x(si);
%
% xi = interp1(xu,x,xiu);
for i=1:length(xiu)
    [dummy ci] = min(abs(xiu(i)-xu));
    xi(i) = x(ci);
end