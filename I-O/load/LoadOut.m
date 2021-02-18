%
% function LoadOut(DBName, FunName, Detail);
%
% Loads ouput structures of 'FunName' function from the list 'DBName'
% into an output structure called Total (includes corresp filenames)
% Saves it into a mat file called:
% 'All.' FunName '.' Detail '.mat'

function LoadOut(DBName, FunName, Detail);

RootDir = pwd;
db = textread(DBName, '%s');

nmydb = length(db);

AccCnt=1;
Files = {};
Total = struct([]);

for ii=1:nmydb
    dn = db{ii};
    fprintf('%s \n', dn);
    cd(dn);

    fname = [dn '.' FunName Detail '.mat'];
    OutArgs =  load(fname);

    Total = CopyStruct(OutArgs,Total,AccCnt);

    Files{AccCnt,1} = dn;
    AccCnt=AccCnt+1;

    cd(RootDir);

end

Total(1).File =[];
for ii=1:nmydb
    Total(ii,1).File = Files(ii);
end


save(['All.' FunName '.' Detail '.mat'], 'Total');