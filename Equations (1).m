

function F = Equations(x)
global Q
  F(1) = x(5)^(2)/(av*x(2)) - ((Q(8))^(2)/Q(2))*exp(-DO2/T) ;
  F(2)=x(4)^(2)/(av*x(1)) - (Q(7)^(2)/Q(1))*exp(-DN2/T) ;
 F(3) =x(5)*x(4)/(av*x(3)) - (Q(8)*Q(7)/Q(3))*exp(-DNO/T) ;
 F(4)= x(11)*x(9)/(av*x(4)) - (Q(11)*Q(9)/Q(7))*exp(-IN/T);
  F(5)=x(11)*x(10)/(av*x(5)) - (Q(11)*Q(10)/Q(8))*exp(-IO/T) ;
 F(6)= x(11)*x(8)/(av*x(3)) - (Q(11)*Q(6)/Q(3))*exp(-INO/T) ;
 F(7)= x(11)*x(7)/(av*x(2)) - (Q(11)*Q(5)/Q(2))*exp(-IO2/T) ;
 F(8)= x(11)*x(6)/(av*x(1)) - (Q(11)*Q(4)/Q(1))*exp(-IN2/T);
  
 F(9) = xnop + 2*xo2p + 2*xn2p + xnp + xop - xe;
 F(10) = 2*xo2 + 2*xo2p + xo + xop + xno + xnop - (2*8.41);
 F(11) = 2*xn2 + 2*xn2p + xn + xnp + xno + xnop - (2*31.67);
 
  

