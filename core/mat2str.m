function  s = mat2str(M,format)

%MAT2STR -- String representation of a numerical array of arbitrary size.
%
%  S = MAT2STR(M)  converts the array M into a sting representation S.
%    S has the form 'reshape([..data..],[size])'  and can be read back
%    by MATLAB to reconstruct the matrix. The default format of each
%    element is '%g'  (see SPRINTF).
%
%  S = MAT2STR(M,FORMAT) uses the format string FORMAT (see SPRINTF).
%
%  Examples:
%    mat2str([1:0.5:3]')  ==>  'reshape([1 1.5 2 2.5 3 ],[5 1 ])'
%    mat2str(magic(3))    ==>  'reshape([8 3 4 1 5 9 6 7 2 ],[3 3 ])'
%    mat2str([])          ==>  'reshape([],[0 0 ])'
%
%  See also MAT2STR1, VECT2STR, NUM2STR, INT2STR, SPRINTF, SSCANF, EVAL.

%  Original coding by Alexander Petrov, UC Irvine.
%  $Revision: 1.0 $  $Date: 2004/03/03 10:30 $

if (nargin==1)  
   format = '%g ' ;
else
   format = [format,' '] ;     % e.g. '%d' becomes '%d '
end

s = ['reshape(',vect2str(M(:)'),',',vect2str(size(M)),')'] ;

%--- Return s
%%%%% End of file MAT2STR.M