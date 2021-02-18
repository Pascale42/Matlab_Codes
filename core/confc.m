function confC = confc(dof,p)
%helper to compute coherence confidence values
if dof <= 2
    confC = 1;
else
    df = 1./((dof/2)-1);
    confC = sqrt(1 - p.^df);
end;
return
