

sets
         f index for bazyaft       /bazyaft1,bazyaft2/
         k index for jodasaz       /joda1,joda2,joda3/
         m index of enhedam        /enhdam1,enhedam2/
         c index of bimarestan     /m22,m9,m6,m15,m2/
         t index of priod         /1,2/
         s                        /s1,s2/;


parameters

         v2(k)           /joda1 150,joda2 250,joda3 500/
         b(t)            /1 0.4,2 0.7/
         ff(t)           /1 700,2 900/
         ttal            /24/
         vv              /0.0001/
         r2              /2/
         mm               /1000000/
         p(s)            /s1 1000,s2 1300/
         RO(T)           /1 2,2 1/
       ;

table d(c,t)

           1       2
m22       30      70

m9        20      50

m6        70      30

m15       20      80

m2        50      50
;

table c1(k,c,t,s)

                         1.s1    1.s2    2.s1    2.s2
joda1.m22                 10      10      10      10

joda1.m9                  10      10      10      10

joda1.m6                  22      22      22      22

joda1.m15                 75      75      75      75

joda1.m2                  55      55      55      55

joda2.m22                 22      22      22      22

joda2.m9                  15      15      15      15

joda2.m6                  10      10      10      10

joda2.m15                 22      22      22      22

joda2.m2                  55      55      55      55

joda3.m22                 30      30      30      30

joda3.m9                  22      22      22      22

joda3.m6                  10      10      10      10

joda3.m15                 22      22      22      22

joda3.m2                  30      30      30      30

;

table c2(m,k,t,S)

                         1.S1    1.S2    2.S1    2.S2
enhdam1.joda1            2       2       2       2

enhdam1.joda2            20      20      20      20

enhdam1.joda3            22      22      22      22

enhedam2.joda1           50      50      50      50

enhedam2.joda2           10      10      10      10

enhedam2.joda3           12      12      12      12
;

table c3(f,k,t,S)

                         1.S1    1.S2    2.S1    2.S2
bazyaft1.joda1           10      10      10      10

bazyaft1.joda2           10      10      10      10

bazyaft1.joda3           15      15      15      15

bazyaft2.joda1           40      40      40      40

bazyaft2.joda2           15      15      15      15

bazyaft2.joda3           15      15      15      15

;

table t1(k,c,S)

             m22.S1         m22.S2          m9.S1          m9.S2          m6.S1              m6.S2            m15.S1       m15.S2               m2.S1       m2.S2
joda1        4              4               4              4               4                 4                  4          4                     4          4

joda2        4              4               3              4               3                 4                  4          4                     4          4

joda3        4              4               3              4               4                 4                  4          4                     4          4

;


table t2(m,k,S)

                 joda1.S1        joda1.S2        joda2.S1        joda2.S2        joda3.S1        joda3.S2

enhdam1            1             1               2               2               3               3

enhedam2           4             4               3               3               4               4
;

table t3(f,k,S)

                 joda1.S1        joda1.S2        joda2.S1        joda2.S2        joda3.S1        joda3.S2
bazyaft1           1             1               1               1               1               1

bazyaft2           2             1               1               1               1               1

;



variables
         q1(k,c,t,S)
         q2(m,k,t,S)
         q3(f,k,t,S)
         a(c,t)
         x1(k,c,t,S)
         x2(m,k,t,S)
         x3(f,k,t,S)
         z
         z1
         z2;
positive variable
PP(T)
QE(C,T)
sss
;

integer variables
         q1(k,c,t,S)
         q2(m,k,t,S)
         q3(f,k,t,S)
         a(c,t)

;

binary variables

x1(k,c,t,S)
x2(m,k,t,S)
x3(f,k,t,S)
;

Equations
         objectiveFunction1
         co1
         CO11
         co2
         co3
         co4
         co5
         co6
         co7
         goal2
         epsilon
;

objectiveFunction1.. sum((c,k,t,S),q1(k,c,t,S)*c1(k,c,t,S))+sum((m,k,t,S),q2(m,k,t,S)*c2(m,k,t,S))+sum((f,k,t,S),q3(f,k,t,S)*c3(f,k,t,S))+sum((c,t),a(c,t)*ff(t))=e=z1;

co1(C,T)..sum((k,S),q1(k,c,t,S))+a(c,t-1)=e=d(c,t);


CO11(C,T)..RO(T)*PP(T)+QE(C,T)=G=d(c,t);

co2(k,t)..sum((c,S),q1(k,c,t,S))=e=sum((m,S),q2(m,k,t,S))/b(t);

co3(k,t)..sum((c,S),q1(k,c,t,S))=e=sum((f,S),q3(f,k,t,S))/(1-b(t));

co4..sum((k,c,t,S),x1(k,c,t,S)*t1(k,c,S))+sum((k,m,t,S),x2(m,k,t,S)*t2(m,k,S))+sum((k,f,t,S),x3(f,k,t,S)*t3(f,k,S))=l=ttal;

co5(k,c,t,S)..q1(k,c,t,S)=l=mm*x1(k,c,t,S)*P(S);

co6(k,m,t,S)..q2(m,k,t,S)=l=mm*x2(m,k,t,S)*P(S);

co7(k,f,t,S)..q3(f,k,t,S)=l=mm*x3(f,k,t,S)*P(S);
goal2..sum((k,c,t,S),x1(k,c,t,S)*t1(k,c,S))+sum((k,m,t,S),x2(m,k,t,S)*t2(m,k,S))+sum((k,f,t,S),x3(f,k,t,S)*t3(f,k,S))=e=z2;
epsilon..z1+z2=e=z;




model supplychain /all/;
option limrow=100, limcol=100;

q1.up(k,c,t,S)=1000;
q2.up(m,k,t,S)=1000;
option MIp=CPlex , optca=0, optcr=0;
solve supplychain minimaizing z using MIP;




















