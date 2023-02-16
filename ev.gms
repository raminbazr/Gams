sets
         i       warehouse                       /i1*i5/
         j       hospital                        /j1*j10/


parameters g(i)     warehouse operating costs       /i1 25000000,i2 20000000,i3 12000000,i4 6000000,i5 12000000/
          l(i)     capacity of warehouse for a medical center    /i1 20000,i2 25000,i3 30000,i4 10000,i5 5000/;

scalar    e     maximum amount available of each medical supply type /500000/;


 variable
         z                 objective
         S(i)              inventory level of medical supply k in warehouse i
         t(i,j)          amount of medical supply k to be delivered from warehouse i to hospital j under disaster scenari
         y(j)            amount of unfulfilled demande
;

positive variable         S,t,y;
binary variable
          x(i)           if warehouse i is selected otherwise 0 ;


equations
         obj           objective
         cons1(i)      constraint for calculating Q
         cons2(j)    constraint for calculating Q
         cons3(j)    constraint for calculating Q
         cons4
         cons5(i)
;

obj..                        z =e= sum(i,g(i)*x(i))+53.49* sum((i,j),t(i,j))+3000* sum((j),y(j));
cons1(i)..                     sum((j),t(i,j))=l=S(i);
cons2(j)..                     sum(i,t(i,j))=e= 4475.595 - y(j);
cons3(j)..                     3000 *y(j)=l=100000000;
cons4..                        sum(i,S(i))=l=e;
cons5(i)..                     S(i)=l=L(i)*x(i);


model location /all/;
solve location using mip minimizing z;
display z.l,x.l,y.l;
