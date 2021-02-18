function Dat = LoadSegments(FileBase, Type, Seg, Channels)
%function Dat = LoadSegments(FileBase, Type, Seg, Channels)

Par= LoadPar([FileBase '.par']);
Dat =[];
nSeg = size(Seg,1);
for i=1:nSeg
    beg = Seg(i,1)*Par.nChannels*2;
    len = diff(Seg(i,:));
    chunk = bload([FileBase '.' Type], [Par.nChannels, len], beg);
    Dat(:,[end+1:end+len]) = chunk(Channels,:);
end
 