%IRS
t0 = 874.5;
n=1170;
Ak(1:694,1)=IRS(1:694,2);
Ak(695:n,1)=IRS(695:n,2)-360;
pAz = polyfit((IRS(1:n,1)-t0),Ak(1:n,1),3);
%RI
%Resat
%Asat
%Msat


pEl = polyfit((IRS(1:n,1)-t0),IRS(1:n,3),3);
pAzr = polyfit((IRS(1:n,1)-t0),IRS(1:n,4),3);
pElr = polyfit((IRS(1:n,1)-t0),IRS(1:n,5),3);

%use to compare graphs
yAz=polyval(pAz,(IRS(1:n,1)-t0));
yEl=polyval(pEl,(IRS(1:n,1)-t0)); 
yAzr=polyval(pAzr,(IRS(1:n,1)-t0)); 
yElr=polyval(pElr,(IRS(1:n,1)-t0)); 

Result(1)=polyval(pAz,0);
Result(2)=polyval(pEl,0);
Result(3)=polyval(pAzr,0);
Result(4)=polyval(pElr,0);

%plot(IRS(1:694,1),IRS(1:694,2),IRS(695:n,1),IRS(695:n,2)-360,IRS(1:n,1),yAz);
%plot(IRS(1:n,1),IRS(1:n,3),IRS(1:n,1),yEl);
%plot(IRS(1:n,1),IRS(1:n,4),IRS(1:n,1),yAzr);
%plot(IRS(1:n,1),IRS(1:n,5),IRS(1:n,1),yElr);

