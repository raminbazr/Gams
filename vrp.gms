$ontext
this model is provide by ramin bazrafshan by student-number 98125030
$offtext
sets
i name of starting city   /1*8/
k name of vehicle      /1*3/
o storage /1/
alias(i,j,p);

table dist(i,j) distance between i & j
          1         2          3          4          5          6           7         8
1         0        990        2160       1060       500        2050        999       2050
2         0        0          1160       1020       590        1060        1220      1050
3         0        1160       0          1740       1720       210         160       250
4         0        1020       1740       0          710        1530        1900      1520
5         0        590        1720       710        0          1580        1820      1570
6         0        1060       210        1530       1580       0           370       30
7         0        1220       160        1900       1820       370         0         400
8         0        1050       250        1520       1570       30          400       0    ;

parameters
q(i) demand of cities / 1 0,2 6,3 3,4  8,5  7,6  9,7  4,8  5/
E(j) earliest time /1 0,2 0,3 0,4  0,5  120,6  120,7  120,8  120/
L(j) latest time /1 0,2 1440,3 1440,4 1440,5 1440,6 1440,7 1440,8 1440/
s(i) time of disembarkation vehicles in cities /1 0,2 20,3 20,4 15,5 20,6 15,7 20,8 15/
;

scalar  aa impact of delay time /0.3/;
scalar vcap capacity of vehicles /18/;
scalar m number of vehicle /3/  ;

variables
TC  ;
binary variable x(i,j,k)
positive variable t(i)
positive variable u(i,k)
positive variable d(j)
positive variable w(i)
;

t.up(i) =1000;
equation
objectivefunction
co1
co2
co3(k)
co4(j)
co5(i)
co6(p,k)
co7(i,j,k)
co8(i,j)
co9(i,j,k)
co10(i,k)
co11(i,j)
co12(j)
co13(j)
;
objectivefunction ..              TC=e= sum((i,j,k), dist(i,j)*x(i,j,k)) + aa*sum (j,d(j)) ;
co1                                    ..sum((i,k),x(i,'1',k))=e=m    ;
co2                                    ..sum((k,j),x('1',j,k))=e=m  ;
co3(k)                                 ..sum(i,q(i)*sum(j,x(i,j,k)))=l= VCAP  ;
co4(j)$(ord(j)<>1)                     ..sum((i,k)$(ord(i)<>ord(j)),x(i,j,k))=e=1   ;
co5(i)$(ord(i)<>1)                     ..sum((j,k)$(ord(i)<>ord(j)),x(i,j,k))=e=1   ;
co6(p,k)                               ..sum(i,x(i,p,k))-sum(j,x(p,j,k))=e=1 ;
co7(i,j,k)$(ord(i)<>1 and ord(j)<>1)   ..(u(i,k)- u(j,k) + VCAP *x(i,j,k)) =l= VCAP - q(j)      ;
co8(i,j)                               ..q(i)+q(j)=l=vcap ;
co9(i,j,k)                             ..q(j)=g= u(i,k)   ;
co10(i,k)                              ..u(i,k)=l= VCAP;
co11(i,j)                              ..(sum((k),x(i,j,k))-1)*1000 +(t(i)+w(i)+ s(i) + dist(i,j)) =l= t(j) ;
co12(j)                                ..E(j) =l= t(j) + w(j) ;
co13(j)                                ..L(j) =g= t(j) - d(j) ;

model vrp /all/       ;
option limrow=100,limcol=100 ;
option  optca=0 , optcr=0;
option mip=cplex;
solve vrp using mip minimizing TC  ;
display TC.l , t.l , u.l, d.l , w.l;



