$ontext
Dantzig-Wolfe Decomposition
$offtext

sets
i 'Reservoir'             /1*3 /
j 'RHS for constraint'   /1*3 /
;

parameter
cost(i)                  /1 20,2 10,3 5/
b(j)                     /1 285,2 360,3 150/
costsub(i)               /1 -20,2 -10,3  5/
;

*-----------------------------------------------------------------------
*                        direct LP formulation
*-----------------------------------------------------------------------

positive variable x(i)   'The amount of water released from reservoir'
free variable z          'objective variable'
;
equations
obj
Demand
limitR4
limitR(i)
;

obj                         .. z =e= sum(i, cost(i)*x(i));
Demand                      .. sum(i, x(i)) =g= b('1');
limitR4                     .. sum(i, x(i)) =l= b('2');
limitR(i)                   .. x(i) =l= b('3');

model power/all/;
solve power minimizing z using lp;

*-----------------------------------------------------------------------
*                                subproblems
*-----------------------------------------------------------------------

positive variables xsub(i);
variables zsub;

parameters
b_sub     'RHS'
c(i)      'cost coefficients'
u1        'dual of Demand'
u2        'dual of limitR4'
v         'dual of convexity constraint'
vp
;

equations
limitR_sub
zsubproblem
;

limitR_sub(i)            .. xsub(i) =l= b_sub;
zsubproblem              .. zsub =e= sum((i,j), (c(i)-(u1+u2))*xsub(i))-vp ;

model sub  /limitR_sub, zsubproblem/;

*-----------------------------------------------------------------------
*                            master problem
*-----------------------------------------------------------------------

set k 'point count' /1*1000/;
set pk(k);
pk(k)=no;
parameter point(i,k);
parameter pointcost(k);
point(i,k) = 0;
pointcost(k) = 0;

positive variables  lambda(k) ;
free variable    zmaster;

equations

obj_master
limitR4_master
Demand_master
convex_master
;

obj_master               .. zmaster =e= sum(pk, pointcost(pk)*lambda(pk));
Demand_master            .. sum((i,pk),point(i,pk)*lambda(pk)) =g= b('1');
limitR4_master           .. sum((i,pk),point(i,pk)*lambda(pk)) =l= b('2');
convex_master            .. sum( pk(k), lambda(k)) =e= 1;

model master  /obj_master, limitR4_master, Demand_master, convex_master/;

*-----------------------------------------------------------------------
*                 options to reduce solver output
*-----------------------------------------------------------------------

option limrow=0;
option limcol=0;

*
*INITIALIZATION PHASE
*

set kk(k) 'current point';
kk('1') = yes;

*
* solve subproblem, check feasibility
*


c(i) = costsub(i) ;
b_sub = b('3');
u1 = 0;
u2 = 0;
v=0;
vp=0;

solve sub using lp minimizing zsub;
point(i,kk) = xsub.l(i);
c(i)=cost(i);
pointcost(kk) = sum((i), c(i)*xsub.l(i));
pk(kk)=yes;
kk(k)=kk(k-1);

display point;
*-----------------------------------------------------------------------
* DANTZIG-WOLFE ALGORITHM
* while (true) do
* solve restricted master
* solve subproblems
* until no more point
*-----------------------------------------------------------------------

set iter 'maximum iterations' /1*10/;
scalars
done /0/
count /0/
phase /1/
iteration
;

*
*iteration
*
loop(iter$(not done),
count = 0;
iteration = ord(iter);

*
* solve master problem to get duals
*

solve master minimizing zmaster using lp;
display zmaster.l,lambda.l,Demand_master.m,limitR4_master.m,convex_master.m;

u1= Demand_master.m;
u2 = limitR4_master.m;
v=  convex_master.m ;

*
* solve each subproblem
*

c(i) = cost(i);
b_sub = b('3');
vp=v;

solve sub using lp minimizing zsub;
display zsub.l,xsub.l;

*
* point
*

if (zsub.l <0 ,
display "new point", count,xsub.l;
point(i,kk) = xsub.l(i);
pointcost(kk) = sum(i, c(i)*xsub.l(i) );
pk(kk)=yes;
kk(k) = kk(k-1);
count = count+ 1;
);
);

*-----------------------------------------------------------------------
*                        recover solution
*-----------------------------------------------------------------------

parameter xfinal(i);
xfinal(i) = sum(pk(k), point(i,pk)*lambda.l(pk));
display xfinal;;
parameter totalcost;
totalcost = sum((i), cost(i)*xfinal(i));
display totalcost;
