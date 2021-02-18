function out = DetetctIntraDir(FileBase)
dn = dir('*.dat');
dati = strfind(dn(1).name,'.dat');
FirstFileName = dn(1).name(1:dati-1);
out = DetectIntraSpk(FirstFileName);

