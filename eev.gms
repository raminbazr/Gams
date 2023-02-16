   sets
         i       warehouse                       /i1*i5/
         j       hospital                        /j1*j10/
         k       senario                         /k1*k3/;

parameters g(i)     warehouse operating costs       /i1 25000000,i2 20000000,i3 12000000,i4 6000000,i5 12000000/
           l(i)     capacity of warehouse for a medical center    /i1 20000,i2 25000,i3 30000,i4 10000,i5 5000/
           pw(k)    probability for each senario       /k1 0.275,k2 0.175,k3 0.55/;




table w(k,j)  cost of location
$CALL =XLS2GMS R=l39:v42 I=E:\industman\twostageexcel.xlsx O=dat1.inc
$include dat1.inc
;
table d(k,j)   capacity of destribution centers
$CALL =XLS2GMS R=l44:v47 I=E:\industman\twostageexcel.xlsx O=dat4.inc
$include dat4.inc
;
table tow(k,j)   capacity of destribution centers
$CALL =XLS2GMS R=l49:v52 I=E:\industman\twostageexcel.xlsx O=dat2.inc
$include dat2.inc
;

scalar    e     maximum amount available of each medical supply type /500000/;

table     c(k,i,j)
$CALL =XLS2GMS R=a39:k54 I=E:\industman\twostageexcel.xlsx O=dat3.inc
$include dat3.inc
;

variable
         z                 objective
         S(i)              inventory level of medical supply k in warehouse i
         t(k,i,j)          amount of medical supply k to be delivered from warehouse i to hospital j under disaster scenari
         y(k,j)            amount of unfulfilled demand
;

positive variable         S,t,y;
binary variable
          x(i)           if warehouse i is selected otherwise 0 ;


equations
         obj           objective
         cons1(i)      constraint for calculating Q
         cons2(k,j)    constraint for calculating Q
         cons3(k,j)    constraint for calculating Q
         cons4
         cons5(i)
         cons6
         cons7
         cons8
         cons9
         cons10

;

obj..                        z =e= sum(i,g(i)*x(i))+sum(k,pw(k))*(sum((k,i,j),c(k,i,j)*t(k,i,j))+ sum((k,j),w(k,j)*y(k,j)));
cons1(i)..                       sum((k,j),t(k,i,j))=l=S(i);
cons2(k,j)..                     sum(i,t(k,i,j))=e= d(k,j)- y(k,j);
cons3(k,j)..                     w(k,j)*y(k,j)=l=tow(k,j);
cons4..                          sum(i,S(i))=l=e;
cons5(i)..                       S(i)=l=L(i)*x(i);
cons6..                          x('i3')=e=1;
cons7..                          x('i4')=e=1;
cons8..                          x('i5')=e=1;
cons9..                          x('i2')=e=0;
cons10..                         x('i1')=e=0;



model location /all/;
solve location using mip minimizing z;
display z.l,y.l;
