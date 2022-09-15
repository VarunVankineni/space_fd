 

lat = (13.0344722)*pi()/180 ;
lon = (77.51169)*pi()/180 ;   %BL1
% lon=(77.51095)*pi()/180 ;   %BL2
T= [julian_date(2016, 11, 16, 18, 29, 0),...
    julian_date(2016, 11, 16, 18, 30, 0),...
    julian_date(2016, 11, 16, 18, 31, 0)] ;
AZ_EL = [273.11 275.53 278.82;6.23 10.53 15.65]*pi()/180; 
% lat = (42+24/60)*pi()/180; 
% lon = -(71+4/60)*pi()/180 ;
% T = 2451545+[0 5 10]/1440 ;
% AZ_EL = [10 70 130; 5 45 10]*pi()/180; 
[rv,L,Ld,Ldd,RVATT] = laplace_orbit_fit(lat,lon,839.91,T,AZ_EL);

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