%problem 3
mu = 398600;
a = 140000;
ec = 0.95;
%specific energy 
E = -mu/(2*a);

%specific angular momentum
h = sqrt(mu*a*(1-ec^2));


%radial velocity, normal velcity and flight path angle 
t = -pi:pi/180:pi;
v_r = mu.*ec.*sin(t)./h;
v_n = mu.*(1+ec.*cos(t))./h;
g = atan((ec.*sin(t))./(1+(ec.*cos(t))));

fprintf('\n %g %g \n',E,h);

%plot(t,v_r,t,v_n);





