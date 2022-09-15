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


tt=Sat(1,:);
%IRS=IRS;
%sidereal time with time base from Az, El measurements
St = LST(tt(1), tt(2), tt(3), IRS(:,1),Elong);%change year month and date here
%right ascension and declination for each instance of time from Az,El
%measurements
[asc, dec]= alphadelta(IRS(:,2),IRS(:,3),St(:,1),lat);
asc=pi+0.5321+(asc*deg);
dec=dec*deg;
%L, Ldot measurements by polynomial fitting of the asc and dec with time
%base in second by taking origin at the mid point of measurement
hue= [julian_date(2016, 11, 16, 17, 29, 0),...
    julian_date(2016, 11, 16, 17, 30, 0),...
    julian_date(2016, 11, 16, 17, 31, 0)] ;
[L,Ld,Ldd] = Lfinder(asc,dec,hue);

%calculate R of the observatory at each St instance
R = [(Re/sqrt(1-(2*f - f*f)*(sind(lat))^2) + H)*cosd(lat)*cosd(St(2)),...
    (Re/sqrt(1-(2*f - f*f)*sind(lat)^2) + H)*cosd(lat)*sind(St(2)),...
    (Re*(1 - f)^2/sqrt(1-(2*f - f*f)*sind(lat)^2) + H)*sind(lat)];
Rd= cross(w,R);
Rdd= cross(w,Rd);
R=R.'*1000;
Rd=Rd.'*1000;
Rdd=Rdd.'*1000;

SMA_E       = 6378137.0;% m
MU_E        = 3.986004418E14; % m^3/s^2
MU_E_ER     = MU_E/(SMA_E*SMA_E*SMA_E); % (ER)^3/s^2
L=L(:,2);
RVA=[R;Rd;Rdd];
RVA = RVA/SMA_E;
Rp = norm(RVA(1:3));
    N = dot(L, RVA(1:3));
    D  = det(horzcat(L, Ld, Ldd))*2;
    D1 = det(horzcat(L, Ld, RVA(7:9)));
    D2 = det(horzcat(L, Ld, RVA(1:3)));
    D3 = det(horzcat(L, RVA(7:9), Ldd));
    D4 = det(horzcat(L, RVA(1:3), Ldd));

    c8 = D*D;
    c6 = ((4*N*D - 4*D1)*D1 - D*D*Rp*Rp);
    c3 = (4*MU_E_ER*D2*(N*D - 2*D1));
    c0 = -4*MU_E_ER*MU_E_ER*D2*D2;

    Z = roots([c8 0 c6 0 0 c3 0 0 c0]);
    for (i = 1:numel(Z))
        if (isreal(Z(i)) && Z(i) > 0)
            r = Z(i);
            break;
        end
    end

    rrr = r*r*r;
    rho = -2*(D1/D) - 2*(MU_E_ER/rrr)*(D2/D);
    rhodot = -(D3/D) - (MU_E_ER/rrr)*(D4/D);

    RV = vertcat(rho*L + RVA(1:3), ...
            rhodot*L + rho*Ld + RVA(4:6))*SMA_E;
    r=[RV(1),RV(2),RV(3)]/1000;
    v=[RV(4),RV(5),RV(6)]/1000;
% 
%      P=dot(L,cross(Ld,Rdd));
%      Q=dot(L,cross(Ld,R));
%      D=dot(L,cross(Ld,Ldd));
%      T=dot(L,cross(Rdd,Ldd));
%      S=dot(L,cross(R,Ldd));
%      P=P/D;
%      Q=mu*Q/D;
%      T=T/(2*D);
%      S=mu*S/(2*D);
%      J=2*dot(L,R);
%      Rsq=dot(R,R);
%      c8=1;
%      c6=-((P^2)-(J*P)+Rsq);
%      c3=-((-J+(2*P))*Q);
%      c0=-(Q^2);
%    
%      
%  Z = roots([c8 0 c6 0 0 c3 0 0 c0]);
%     for (i = 1:numel(Z))
%         if (isreal(Z(i)) && Z(i) > 0)
%             res = Z(i);
%             break;
%         end
%     end  
%      
% ro=-(P+(Q/(res^3)));
% rod=-(T+(S/(res^3)));
% 
% rf=(R+ro);
% Ld=Ld.';
% v=Rd+(ro*Ld)+(rod*L);
 vp=norm(v);
 rfp=norm(r);
% Rp=norm(R);






%...Input data:
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
    T2 = 2*pi/sqrt(mu)*coe(7)^1.5; % Equation 2.73
    fprintf('\n Period:')
    fprintf('\n Seconds = %g', T2)
    fprintf('\n Minutes = %g', T2/60)
    fprintf('\n Hours = %g', T2/3600)
    fprintf('\n Days = %g', T2/24/3600)
end
fprintf('\n-----------------------------------------------\n')
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜







