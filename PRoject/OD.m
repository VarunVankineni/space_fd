function [  ] = OD(n)
% wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
%{
deg - factor for converting between degrees and radians
pi - 3.1415926...
mu - gravitational parameter (km^3/s^2)
r1, r2, r3 - three coplanar geocentric position vectors (km)
ierr - 0 if r1, r2, r3 are found to be coplanar
1 otherwise
v2 - the velocity corresponding to r2 (km/s)
coe - orbital elements [h e RA incl w TA a]
where h = angular momentum (km^2/s)
e = eccentricity
RA = right ascension of the ascending node (rad)
incl = orbit inclination (rad)
w = argument of perigee (rad)
TA = true anomaly (rad)
a = semimajor axis (km)
T - period of elliptic orbit (s)
%}

deg = pi/180;
global mu
mu = 398600;

prompt='Choose tracking station(1 or 2):  ';
m=input(prompt);
switch (m)
case 1
    a=[.83991 
        360-282.48831 
        13.0344772];
case 2
        otherwise

    a=[0.83991 
        360-282.48905 
        13.0344722];
end
prompt='Enter time(y, m, d, ut(in hrs)): ';
t=input(prompt);
    a(2)=LST(t(1),t(2),t(3),t(4),a(2));

switch (n)
    case 1
    %%Gibbs method
   
    prompt='Enter three coplanar position vectors in SEZ frame (rho,Az,El):\nr1 = ';
    r1=input(prompt);
    r1=convert(r1(1),r1(2),r1(3),a(1),a(2),a(3));
    prompt='r2 = ';
    r2=input(prompt);
    r2=convert(r2(1),r2(2),r2(3),a(1),a(2),a(3));
    prompt='r3 = ';
    r3=input(prompt);
    r3=convert(r3(1),r3(2),r3(3),a(1),a(2),a(3));
    [v2, ierr] = gibbs(r1, r2, r3);
    %...
%...Echo the input data to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n\n Input data:\n')
fprintf('\n Gravitational parameter (km^3/s^2) = %g\n', mu)
fprintf('\n r1 (km) = [%g %g %g]', r1(1), r1(2), r1(3))
fprintf('\n r2 (km) = [%g %g %g]', r2(1), r2(2), r2(3))
fprintf('\n r3 (km) = [%g %g %g]', r3(1), r3(2), r3(3))    
fprintf('\n\n');


%...If the vectors r1, r2, r3, are not coplanar, abort:
if ierr == 1 && n==1
fprintf('\n These vectors are not coplanar.\n\n')
return
end
coe = coe_from_sv(r2,v2,mu);
h = coe(1);
e = coe(2);
RA = coe(3);
incl = coe(4);
w = coe(5);
TA = coe(6);
a = coe(7);
%...Output the results to the command window:
fprintf(' Solution:')
fprintf('\n');
fprintf('\n v2 (km/s) = [%g %g %g]', v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:');
fprintf('\n Angular momentum (km^2/s) = %g', h)
fprintf('\n Eccentricity = %g', e)
fprintf('\n Inclination (deg) = %g', incl/deg)
fprintf('\n RA of ascending node (deg) = %g', RA/deg)
fprintf('\n Argument of perigee (deg) = %g', w/deg)
fprintf('\n True anomaly (deg) = %g', TA/deg)
fprintf('\n Semimajor axis (km) = %g', a)
%...If the orbit is an ellipse, output the period:
if e < 1
T = 2*pi/sqrt(mu)*coe(7)^1.5;
fprintf('\n Period (s) = %g', T)
end
fprintf('\n-----------------------------------------------------\n')

    case 2
    %%Lamberts method
    prompt='Enter two position vectors in SEZ frame(rho,Az,El), time interval(in s) and type of orbit:\nr1 = ';
    r1=input(prompt);
    r1=convert(r1(1),r1(2),r1(3),a(1),a(2),a(3));
    prompt='r2 = ';
    r2=input(prompt);
    r2=convert(r2(1),r2(2),r2(3),a(1),a(2),a(3));
    prompt='dt = ';
    dt=input(prompt);
    prompt='Type of orbit(pro/retro)= ';
    string=input(prompt);
    [v1, v2] = lambert(r1, r2, dt, string);
    coe = coe_from_sv(r1, v1, mu);
TA1 = coe(6);
coe = coe_from_sv(r2, v2, mu);
TA2 = coe(6);
%...Echo the input data and output the results to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Example 5.2: Lambert”s Problem\n')
fprintf('\n\n Input data:\n');
fprintf('\n Gravitational parameter (km^3/s^2) = %g\n', mu);
fprintf('\n r1 (km) = [%g %g %g]', ...
r1(1), r1(2), r1(3))
fprintf('\n r2 (km) = [%g %g %g]', ...
r2(1), r2(2), r2(3))
fprintf('\n Elapsed time (s) = %g', dt);
fprintf('\n\n Solution:\n')
fprintf('\n v1 (km/s) = [%g %g %g]', ...
v1(1), v1(2), v1(3))
fprintf('\n v2 (km/s) = [%g %g %g]', ...
v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:')
fprintf('\n Angular momentum (km^2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n Inclination (deg) = %g', coe(4)/deg)
fprintf('\n RA of ascending node (deg) = %g', coe(3)/deg)
fprintf('\n Argument of perigee (deg) = %g', coe(5)/deg)
fprintf('\n True anomaly initial (deg) = %g', TA1/deg)
fprintf('\n True anomaly final (deg) = %g', TA2/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))
fprintf('\n Periapse radius (km) = %g', coe(1)^2/mu/(1 + coe(2)))
%...If the orbit is an ellipse, output its period:
if coe(2)<1
T = 2*pi/sqrt(mu)*coe(7)^1.5;
fprintf('\n Period:')
fprintf('\n Seconds = %g', T)
fprintf('\n Minutes = %g', T/60)
fprintf('\n Hours = %g', T/3600)
fprintf('\n Days = %g', T/24/3600)
end
fprintf('\n-----------------------------------------------------\n')
% wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
    
    
    case 3
    case 4
    otherwise
    fprintf('\n Please enter a valid number n=[1,4].\n\n')
   
end
end