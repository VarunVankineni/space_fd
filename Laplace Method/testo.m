
SAT=IRS;
date=SAT(1,6:8);
lat = (13.0344722)*pi()/180 ;
lon = (77.51169)*pi()/180 ;   %BL1
%lon=(77.51095)*pi()/180 ;   %BL2
T= [julian_date(date(1), date(2), date(3), SAT(1,4), SAT(1,5), 0),...
    julian_date(date(1), date(2), date(3), SAT(2,4), SAT(2,5), 0),...
    julian_date(date(1), date(2), date(3), SAT(3,4), SAT(3,5), 0)] ;
AZ_EL = [SAT(1,2) SAT(2,2) SAT(3,2);SAT(1,3) SAT(2,3) SAT(3,3)]*pi()/180;
alt=839.91;
  
constants;
    for (i = 1:numel(T))
        dv = datevec(T(i));
        [asc(1,i),dec(1,i)] = horizontal_to_equatorial(AZ_EL(1,i), ...
                AZ_EL(2,i), lat, lon, dv(1), dv(2), dv(3), dv(4), dv(5), dv(6));
    end

    L = [cos(asc(1,:)).*cos(dec(1,:));
         sin(asc(1,:)).*cos(dec(1,:)); ...
         sin(dec(1,:))];

    Ld = ((T(2)-T(3)) / ((T(1)-T(2))*(T(1)-T(3))) * L(:,1) + ...
            (2*T(2)-T(1)-T(3)) / ((T(2)-T(1))*(T(2)-T(3))) * L(:,2) + ...
            (T(2)-T(1)) / ((T(3)-T(1))*(T(3)-T(2))) * L(:,3))/86400;

    Ldd =(2 / ((T(1)-T(2))*(T(1)-T(3))) * L(:,1) + ...
            2 / ((T(2)-T(1))*(T(2)-T(3))) * L(:,2) + ...
            2 / ((T(3)-T(1))*(T(3)-T(2))) * L(:,3))/(86400*86400); 
        
    dv = datevec(T(2));
    [RVA,~] = geodetic_to_ECI(lat, lon, alt, dv(1), dv(2), dv(3), dv(4), dv(5), dv(6));
    RVA = RVA/SMA_E;

    R = norm(RVA(1:3));
    N = dot(L(:,2), RVA(1:3));
    D  = det(horzcat(L(:,2), Ld, Ldd))*2;
    D1 = det(horzcat(L(:,2), Ld, RVA(7:9)));
    D2 = det(horzcat(L(:,2), Ld, RVA(1:3)));
    D3 = det(horzcat(L(:,2), RVA(7:9), Ldd));
    D4 = det(horzcat(L(:,2), RVA(1:3), Ldd));

    c8 = D*D;
    c6 = ((4*N*D - 4*D1)*D1 - D*D*R*R);
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

    RV = vertcat(rho*L(:,2) + RVA(1:3), ...
            rhodot*L(:,2) + rho*Ld + RVA(4:6))*SMA_E;
rv=RV;
global mu;
deg = pi/180;
mu = 398600;

%...Input data:
r = [rv(1) rv(2) rv(3)]/1000;
rp=norm(r);
v = [rv(4) rv(5) rv(6)]/1000;
vp=norm(v);
%...
%...Algorithm 4.1:
coe = coe_from_sv(r,v);
%...Echo the input data and output results to the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 4.3\n')
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
        
        
        
        
        