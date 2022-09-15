%constants
mu = 398600;%km3/sec2
J2 = 1.082627e-3;
Re=6378.145;

%params
a=42164.2;
i=3*pi()/180;
e=0.04;

n=sqrt(mu/a^3);%frequency
p=a*(1-e^2);%semi latusrectum
RANs= -1.5*n*J2*cos(i)*((Re/p)^2);%rad/sec
AoPs= 0.75*n*J2*((5*(cos(i)^2))-1)*((Re/p)^2);%rad/sec
RANd= RANs*86400*180/(pi());%degrees/day
AoPd= AoPs*86400*180/(pi());%degrees/day