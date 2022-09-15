function [L, Ld,Ldt,yasc,ydec,yascd,ydecd,pasc,pascdot] = Lunits(t,asc,dec)
t0=(range(t)/2)+t(1);
pasc = polyfit((t-t0),asc,4);
pascdot=polyder(pasc);
yasc = polyval(pasc,t-t0);
yascd = (72.92e-6)+(polyval(pascdot,t-t0)/(180/pi));
pdec = polyfit((t-t0),dec,4);
pdecdot = polyder(pdec);
ydec = polyval(pdec,t-t0);
ydecd = polyval(pdecdot,t-t0)/(180/pi);
L = [(cosd(yasc).*cosd(ydec)) ,(sind(yasc).*cosd(ydec)), (sind(ydec))];
Ld= [(-(sind(yasc).*cosd(ydec).*yascd)-(cosd(yasc).*sind(ydec).*ydecd)), ((cosd(yasc).*cosd(ydec).*yascd)-(sind(yasc).*sind(ydec).*ydecd)), (cosd(ydec).*ydecd)];
 [m,n]=size(L);
%  for i=2:(m-1)
%      Ld(i,:)=(L(i+1,:)-L(i-1,:))/(t(i+1)-t(i-1));
%  end
%  Ld(1,:)=Ld(2,:);
%  Ld(m,:)=Ld(m-1,:);
 for i=2:(m-1)
     Ldt(i,:)=(Ld(i+1,:)-Ld(i-1,:))/(t(i+1)-t(i-1));
 end
 Ldt(1,:)=Ldt(2,:);
 Ldt(m,:)=Ldt(m-1,:);




