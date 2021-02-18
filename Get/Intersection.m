%  function [C,IA, IB] = Intersection(A, B)
%
% find the intersection (C) between rows of matrices A and B
% as well as their index position in A (IA) and B (IB)
%
% A and B are matrices of n rows x 2 columns


function [C,IA, IB] = Intersection(A, B)

C=[];
IA=[];
IB=[];

if size(A) == [1 2]
lengthA=1;
else
    lengthA=length(A);
end
if size(B)== [1 2]
    lengthB=1;
else
lengthB=length(B);
end

for i=1:lengthA
    for j=1:lengthB
        if A(i,1)==B(j,1) & A(i,2)==B(j,2)
            c=B(j,:); C=[C;c];
            ia=find(A(:,1) == A(i,1) & A(:,2) == A(i,2));
            IA=[IA;ia];
            ib=find(B(:,1) == B(j,1) & B(:,2) == B(j,2));
            IB=[IB;ib];
        end
    end
end
