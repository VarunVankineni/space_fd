%constants
global mu
deg = pi/180;
mu = 398600;
f = 1/298.256421867;
Re = 6378.13655;
H= 0.83991;
w= [0, 0, 72.9217e-6];
lat=13.0344722;
Elong=(360- 282.48905);
%Input ro and rodot for the mid point to compute r and rdot 
ro=907.1252103;
rod=-3.71756712;
[k1,k2] = size(IRS);
%k=fix(k1/2);
k=340;
%sidereal time with time base from Az, El measurements
St = LST(2016, 9, 8, IRS(:,1),Elong);%change year month and date here
%right ascension and declination for each instance of time from Az,El
%measurements
[asc, dec]= alphadelta(IRS(:,2),IRS(:,3),St(:,1),lat);
%L, Ldot measurements by polynomial fitting of the asc and dec with time
%base in second by taking origin at the mid point of measurement
[L,Ld,Ldt,yasc,ydec,yascd,ydecd,pasc,pasdot]= Lunits(IRS(:,6),asc,dec);

%calculate R of the observatory at each St instance
for i = 1:k1
    R(i,:) = [(Re/sqrt(1-(2*f - f*f)*(sind(lat))^2) + H)*cosd(lat)*cosd(St(i)),(Re/sqrt(1-(2*f - f*f)*sind(lat)^2) + H)*cosd(lat)*sind(St(i)),(Re*(1 - f)^2/sqrt(1-(2*f - f*f)*sind(lat)^2) + H)*sind(lat)];
    Rd(i,:)= cross(w,R(i,:));
    Rdd(i,:)= cross(w,Rd(i,:));
end
for i=1:(k1-1)
          Rd(i,:)=(R(i+1,:)-R(i,:))/(IRS(i+1,6)-IRS(i,6));
end
      Rd(k1,:)=Rd(k1-1,:);
for i=1:(k1-1)
          Rdd(i,:)=(Rd(i+1,:)-Rd(i,:))/(IRS(i+1,6)-IRS(i,6));
end
      Rdd(k1,:)=Rdd(k1-1,:);
 for i = 1:k1
     P(i,1)=dot(L(i,:),cross(Ld(i,:),Rdd(i,:)));
     Q(i,1)=dot(L(i,:),cross(Ld(i,:),R(i,:)));
     D(i,1)=dot(L(i,:),cross(Ld(i,:),Ldt(i,:)));
     T(i,1)=dot(L(i,:),cross(Rdd(i,:),Ldt(i,:)));
     S(i,1)=dot(L(i,:),cross(R(i,:),Ldt(i,:)));
     P(i,1)=P(i,1)/D(i,1);
     Q(i,1)=mu*Q(i,1)/D(i,1);
     T(i,1)=T(i,1)/(2*D(i,1));
     S(i,1)=mu*S(i,1)/(2*D(i,1));
     J(i,1)=2*dot(L(i,:),R(i,:));
     Rsq(i,1)=dot(R(i,:),R(i,:));
     c8(i,1)=1;
     c6(i,1)=-((P(i)^2)+(J(i)*P(i))+Rsq(i));
     c3(i,1)=-((J(i)+(2*P(i)))*Q(i));
     c0(i,1)=-(Q(i)^2);
 end
 for i=1:k1-2
    Z(i,:) = roots([c8(i) 0 c6(i) 0 0 c3(i) 0 0 c0(i)]);
     for j = 1:numel(Z(i,:))
         if (isreal(Z(i,j)) && Z(i,j) > 0)
             res(i,1) = Z(i,j);
             break;
         end
     end
     
 end
 %calculate r for the mid point
pika=;
 ro=P(pika)+(Q(pika)/(res(pika)^3));
 rod=T(pika)+(S(pika)/(res(pika)^3));
rf=(R(k,:)+ ro*L(k,:));
v=Rd(k,:)+(ro*Ld(k,:))+(rod*L(k,:));
vp=norm(v);
rfp=norm(rf);
Rk=(R(k,:));
Rp=norm(R(k,:));
%placeholders for plotting to check accuracy
%plot(IRS(1:1170,1),asc(:,1),IRS(1:1170,1),yasc(:,1));
%plot(IRS(1:1170,1),asc,IRS(1:1170,1),yasc);


%...Input data:
r = rf;
%...
%...Algorithm 4.1:
coe = coe_from_sv(r,v);
%...Echo the input data and output results to the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 4.3\n')
fprintf('\n Mid point = %g', k)
%fprintf('\n Gravitational parameter (kmˆ3/sˆ2) = %g\n', mu)
fprintf('\n State vector:\n')
fprintf('\n r (km) = [%g %g %g]', ...
r(1), r(2), r(3))
fprintf('\n v (km/s) = [%g %g %g]', ...
v(1), v(2), v(3))
disp(' ')
fprintf('\n Angular momentum (kmˆ2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n Right ascension (deg) = %g', coe(3)/deg)
fprintf('\n Inclination (deg) = %g', coe(4)/deg)
fprintf('\n Argument of perigee (deg) = %g', coe(5)/deg)
fprintf('\n True anomaly (deg) = %g', coe(6)/deg)
fprintf('\n Semimajor axis (km): = %g', coe(7))
%...if the orbit is an ellipse, output its period:
if coe(2)<1
    T = 2*pi/sqrt(mu)*coe(7)^1.5; % Equation 2.73
    fprintf('\n Period:')
    fprintf('\n Seconds = %g', T)
    fprintf('\n Minutes = %g', T/60)
    fprintf('\n Hours = %g', T/3600)
    fprintf('\n Days = %g', T/24/3600)
end
fprintf('\n-----------------------------------------------\n')
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜







