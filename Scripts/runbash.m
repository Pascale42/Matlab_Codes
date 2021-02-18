function runbash(DbName)

db = LoadStringArray(DbName);

nmydb = length(db);

RootDir = pwd;

for i=1:nmydb
    dn = db{i}; %current database entry - name of directory
    fprintf('%s \n', dn);
    cd(dn);
    PhaseModulation(db(i), 'compute', 01, 'THE',[3 7], 1);
    cd(RootDir);
end