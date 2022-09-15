function [ alfa ,delta ] = R2RA_Dec( R )
r = norm(R);
l  = R(1)/r; m = R(2)/r; n =R(3)/r;        % Direction cosines
delta = asin(n)*180/pi;                    % Declination
% Right ascension:
if (m >0)
    alfa = acos(l/cosd(delta))*180/pi;
else
    alfa = 360 - acos(l/cosd(delta))*180/pi;
end



