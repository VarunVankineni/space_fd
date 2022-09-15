function[asc, dec]= alphadelta(Az,El,St,lat)
dec(:,1)=asind((cosd(lat)*cosd(Az).*cosd(El))+(sind(lat)*sind(El)));
[m,n] = size(Az);
for i = 1:m
    if (Az(i)<=180)&(Az(i)>=0)
        asc(i,1)=St(i)+acosd(((cosd(lat)*sind(El(i)))-(sind(lat)*cosd(El(i))*cosd(Az(i))))/(cosd(dec(i))));
    elseif(Az(i)>=180)&(Az(i)<=360)
        asc(i,1)=St(i)-acosd(((cosd(lat)*sind(El(i)))-(sind(lat)*cosd(El(i))*cosd(Az(i))))/(cosd(dec(i))));
    end
end

