%%
%Defining the constants and evaluating the quantities which are independent
%of temp needed for finding the partition functions.
global Q
%Characteristic Temperatures matrix of diatomics.
D = zeros(6,2);      
D(1,1) = 2*2.9; D(2,1) = 2*2.1; D(3,1) = 2.44; D(4,1) = D(1,1); 
D(5,1) = D(2,1); D(6,1) = 2.86; %Rotational temperatures
D(1,2) = 3390; D(2,2) = 2270; D(3,2) = 2740; D(4,2) = D(1,2);
D(5,2) = D(2,2); D(6,2) = 3419; %Vibrational temperatures

const = ((2*pi*1.38*10^(-23))^(1.5))*10^(34*3)/(6.626^(3));
q = zeros(11,1);
q(1,1) = (28*1.66*10^(-27))^(1.5);
q(2,1) = (32*1.66*10^(-27))^(1.5);
q(3,1) = (30*1.66*10^(-27))^(1.5);
q(4,1) = q(1);
q(5,1) = q(2);
q(6,1) = q(3);
q(7,1) = (14*1.66*10^(-27))^(1.5);
q(8,1) = (16*1.66*10^(-27))^(1.5);
q(9,1) = q(7);
q(10,1) = q(8);
q(11,1) = (9.1*10^(-31))^(1.5);

QM = zeros(11,4);      %Matrix of different partition functions for each
                       %species.

for i = 7:11          %Monoatomic species wont have rot and vib.
    QM(i,1) = 1;
    QM(i,2) = 1;
end
%%

for i= 1:770
%Doing it at just one temperature at the moment
T = (i*10)+300;
byT = 1/T;

for i = 1:6
    QM(i,1) = T/D(i,1);
    QM(i,2) = 1/(1-exp(-D(i,2)*byT));
end

%Electronic partition functions
QM(1,3) = 1 + 2*exp(-99600*byT);
QM(2,3) = 3 + 2*exp(-11400*byT);
QM(3,3) = 2 + 2*exp(-174*byT) + 2*exp(-63300*byT);
QM(4,3) = QM(1,3);
QM(5,3) = QM(2,3);
QM(6,3) = 2 + 2*exp(-174*byT) + 2*exp(-75090*byT);
QM(7,3) = 4 + 6*exp(-27658.7*byT) + 4*exp(-27671*byT) + 6*exp(-41492.4*byT);
QM(8,3) = 5 + 3*exp(-228*byT) + exp(-326*byT) + 5*exp(-22850*byT);
QM(9,3) = 9 + 5*exp(-22055*byT) + exp(-47070*byT);
QM(10,3) = 4 + 10*exp(-38620*byT) + 6*exp(-58270*byT);
QM(11,3) = 2;

%Evaluating translational partition functions
for i = 1:11
    QM(i,4) = const*q(i)*T^(1.5);
end

%Overall partition functions for all species.
for i = 1:11
    Q(i) = QM(i,1)*QM(i,2)*QM(i,3)*QM(i,4);
end

%We have the Q's for each species.

%%
  DO2 =59500 ;DN2 =113000 ;DNO = 75500;
  IO2 = 142000;IN2 = 181000;INO =108000 ;IN = 169000;IO = 158000;
  av = 6.023*10^(23);
end  
 
  
  






                                        
                                        

