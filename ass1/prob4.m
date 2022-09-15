%problem 4
mu = 398600;
r = 402000;
v = 2.23;
tr = 150*pi/180; %converting to radians


%semi major
a = mu/(v^2 - (2*mu/r));

%solveing for eccentricity
p = [a r*cos(tr) -r-a];
ec=roots(p);

fprintf('\n eccentricity values %g  \n',ec);

%choose e accordingly , here its a value greater than 1 for hyperbola

%closest approach and velocity at the same
r_cl = a*(ec-1);
v_cl = sqrt(mu*((2/r_cl)+(1/a)));

fprintf('\n distance of closest approach and corresponding velocity are %g km %g km/s\n',r_cl(1),v_cl(1));
%use 1 or 2 accordingly

%v_excess
v_x = sqrt(mu/a);
fprintf('\n excess velocity is %g km/s\n',v_x);

%hyperbola
xh = -402000:5:402000;
yh = sqrt(a^2+(xh.^2./(ec(1)^2-1)));
plot(xh,yh);




