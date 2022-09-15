%Solver for the composition at 2500K.(Questions 1a,c,e)
global Q
kp1 = 2.917*10^(-7);kp2 = 0.014487;kp3 = 0.05929;
% = @root2de;%change the 'e' in 2de to the letter of the
Q = 0.05929; 

fun = @root2de;                 %question.
x0 = [1,1];                   %Initial solution estimate.
x = fsolve(fun,x0);

pN = kp1*(x(1))^(0.5);
pO = kp2*(x(2))^(0.5);
pNO = kp3*(x(1)*x(2))^(0.5);

x(1)+x(2)+pO+pNO+ pN