function  s = vect2str(V,format)

%VECT2STR -- String representation of a numerical vector.
%
%  S = VECT2STR(V)  converts the vector V into a sting representation S.
%    S is enclosed in square brackets and can be read back by MATLAB to
%    construct an equivalent vector. Matrices are treated as V(:).
%    The default format of each element is '%g'  (see SPRINTF).
%
%  S = VECT2STR(V,FORMAT) uses the format string FORMAT (see SPRINTF).
%
%  Examples:
%    vect2str([1:0.5:3])   ==>  '[1 1.5 2 2.5 3 ]'
%    vect2str([1:0.5:3]')  ==>  '[1 1.5 2 2.5 3 ]''
%    vect2str(magic(3))    ==>  '[8 3 4 1 5 9 6 7 2 ]''
%    vect2str([])          ==>  '[]'
%
%  See also MAT2STR, NUM2STR, INT2STR, SPRINTF, SSCANF, EVAL.

%  Original coding by Alexander Petrov, UC Irvine.
%  $Revision: 1.1 $  $Date: 2005/03/11 10:20 $

if (nargin==1)  
   format = '%g ' ;
else
   format = [format,' '] ;     % e.g. '%d' becomes '%d '
end

if (size(V,1)<=1)   % row vector or scalar
   s = ['[',sprintf(format,V(:)),']'] ;
else   % column vector (or matrix)
   s = ['[',sprintf(format,V(:)),']'''] ;
end   

%--- Return s
%%%%% End of file VECT2STR.M