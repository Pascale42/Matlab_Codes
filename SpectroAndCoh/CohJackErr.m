%function [confC,phistd,Cerr,C] = CohJackErr(Cxy, Phxy, pval)  
function [confC,phistd,Cerr,C] = CohJackErr(Cxy, Phxy, pval)  
K=1; % for now - for average cross-spectrum acrros tapers.
nf=size(Cxy,2);
nt = size(Cxy,1);
nCh = size(Cxy,3);
pp=1-pval/2;

dim=K*nt;
dof=2*dim;
  
if dof <= 2
   confC = 1;
else     
   df = 1./((dof/2)-1);
   confC = sqrt(1 - pval.^df);
end;
tcrit=tinv(pp,dof-1);


for ch1=1:nCh-1
    for ch2=ch1+1:nCh
        Sxx = sq(Cxy(:,:,ch1,ch1));
        Syy = sq(Cxy(:,:,ch2,ch2));
        Sxy = sq(Cxy(:,:,ch1,ch2)) .*sqrt(Sxx.*Syy);
        Sxy = Sxy.*exp(sqrt(-1)*sq(Phxy(:,:,ch1,ch2)));
        C(:,ch1,ch2)  =abs( sum(Sxy)./sqrt(sum(Sxx).*sum(Syy)));
        for k=1:dim;
            indxk=setdiff(1:dim,k);
         
            mSxx = sq(sum(Sxx(indxk,:)));
            mSyy = sq(sum(Syy(indxk,:)));
            Cxyk=sum(Sxy(indxk,:))./sqrt(mSxx.*mSyy);
            absCxyk=abs(Cxyk);
            atanhCxyk(k,:)=sqrt(2*dim-2)*atanh(absCxyk);
            phasefactorxyk(k,:)=Cxyk./absCxyk;

        end;
        
        atanhC=sqrt(2*dim-2)*atanh(C(:,ch1,ch2));
        sigma12=sqrt(dim-1)*squeeze(std(atanhCxyk,1,1));
        Cu=atanhC+tcrit.*sigma12';
        Cl=atanhC-tcrit.*sigma12';
        Cerr(:,ch1,ch2,1) = tanh(Cl/sqrt(2*dim-2));
        Cerr(:,ch1,ch2,2) = tanh(Cu/sqrt(2*dim-2));
        phistd(:,ch1,ch2)=(2*dim-2)*(1-abs(squeeze(mean(phasefactorxyk))'));

    end;
end
for ch1=1:nCh
  tcrit=tinv(pp,dim-1);
   Sxx = sq(Cxy(:,:,ch1,ch1));
   C(:,ch1,ch1) = (sum(Sxx)/(dim-1))';
   for k=1:dim;
       indxk=setdiff(1:dim,k);
      
       Sjk(k,:)=sq(sum(Sxx(indxk,:)))/(dim-1); % 1-drop spectrum
   end;
   sigma=sqrt(dim-1)*squeeze(std(log(Sjk),1,1)); if C==1; sigma=sigma'; end; 
   conf=tcrit.*sigma;
   conf=squeeze(conf); 
   Cerr(:,ch1,ch1,1)=C(:,ch1,ch1).*exp(-conf'); 
   Cerr(:,ch1,ch1,2)=C(:,ch1,ch1).*exp(conf');
end