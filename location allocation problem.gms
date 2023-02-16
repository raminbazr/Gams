* Benders Decomposition algorithm
set i /i1*i30/;
set j /j1*j20/;
set iter/iter1*iter92/;

parameter  D(i) 'the customers demand'
/
i1  54
i2  53
i3  78
i4  36
i5  40
i6  36
i7  63
i8  61
i9  42
i10 54
i11 49
i12 40
i13 31
i14 75
i15 51
i16 37
i17 77
i18 37
i19 74
i20 35
/;

table c(j,i) 'transportation cost'
           i1        i2         i3         i4         i5         i6         i7         i8         i9         i10        i11        i12        i13        i14        i15        i16        i17        i18        i19        i20
j1        140        531        304        303        542        431        750        443        451        824        186        774        612        458        533        572        243        334        835        663
j2        682        112        559        505        68         201        177        322        270        244        406        234        59         461        345        196        369        415        243        499
j3        729        379        708        660        247        217        439        161        501        480        448        259        247        682        93         479        366        589        401        773
j4        697        120        519        467        222        331        149        470        212        216        460        362        222        372        506        80         453        376        294        352
j5        910        331        787        733        272        391        182        462        488        171        628        104        201        672        437        381        578        642        54         661
j6        373        444        455        427        398        266        638        223        439        706        173        593        461        537        306        520        118        405        686        713
j7        520        371        270        245        497        521        513        652        227        573        443        712        540        82         727        297        502        213        661        129
j8        294        310        205        162        358        291        530        379        212        605        108        606        427        258        468        336        181        125        635        446
j9        740        285        515        473        408        500        293        643        268        332        561        542        414        327        689        186        581        397        444        205
j10       681        395        433        404        528        589        465        730        300        509        571        701        552        240        792        297        616        357        617        38
j11       341        372        102        64         469        445        569        554        216        639        292        712        529        100        640        343        367        73         701        308
j12       750        211        549        500        316        425        182        564        258        227        533        430        311        379        600        134        536        414        333        305
j13       547        65         418        365        128        168        281        311        148        355        285        370        186        338        369        149        269        276        381        419
j14       634        388        642        599        272        182        500        60         479        552        364        361        300        646        32         489        277        538        490        765
j15       467        136        333        280        208        202        357        336        99         431        224        454        271        270        408        177        233        192        465        388
j16       374        296        359        319        281        170        504        225        282        576        86         511        350        397        316        365        41         272        575        560
j17       77         538        244        253        566        467        759        498        438        833        213        806        637        414        589        566        287        301        855        622
j18       845        242        681        628        265        398        22         514        373         54        587        265        218        537        521        244        560        537        145        493
j19       549        247        512        465        156        31         400        112        327        464        262        345        212        497        177        345        189        395        434        613
j20       739        248        659        607        124        198        268        256        397        312        451        144        87         589        242        340        388        522        247        643      

;

parameter f(j)  'fixed cost of facilities'
/
j1  366
j2  314
j3  504
j4  504
j5  580
j6  415
j7  456
j8  549
j9  310
j10 316
j11 459
j12 501
j13 302
j14 415
j15 320
j16 425
j17 506
j18 477
j19 579
j20 554
/;

variables
obj  'the total cost'
x(i,j)  'allocation variables'
y(j) 'location variables'

positive variable x
binary variable y

equations
cost    'the obective function'
eq1(i)   'customers allocation to facilities'
eq2(i,j) 'facilities location constraint' ;

cost .. obj=e=sum((i,j),c(j,i)*D(i)*x(i,j))+sum(j,y(j)*f(j));
eq1(i) .. sum(j,x(i,j))=g=1;
eq2(i,j) .. x(i,j)=l=y(j);

model originalproblem /cost,eq1,eq2/
option optcr=0;
option mip=cplex;
solve originalproblem minimizing obj using mip;
display obj.l, y.l, x.l;

*/////////////////////// Benders Dual Sub Problem ////////////////
y.l(j)=1;
variables
v(i)  'dual variable for allocation constraint'
w(i,j) 'dual variable for location constraint'
zdsp 'the objective function of DSP'

positive variable v, w;

equations
dspobj  'the equation of DSP objective function'
eq3(i,j) 'the DSP constraint'
;

dspobj .. zdsp=e=sum(i,v(i))-sum((i,j),w(i,j)*y.l(j))+sum(j,y.l(j)*f(j));
eq3(i,j) .. v(i)-w(i,j)=l=c(j,i)*D(i);

model dualsubproblem /dspobj, eq3/;

*///////////////////// Benders dual sub problem for unbounded situtions ///////////

equations
extremerayobj      'modified objective function'
extremeraycons     ' modified constraint' ;

extremerayobj .. sum(i,v(i))-sum((i,j),w(i,j)*y.l(j))=e=1;
extremeraycons(i,j) ..  v(i)-w(i,j)=l=0;

model extremeray /extremerayobj,extremeraycons,dspobj/;

*//////////////////// Benders main problem /////////////////////
set optcut(iter)  'dynamic optimal cut set'
set unbcut(iter) 'dynamic feasibility cut set'
;
optcut(iter)=no;
unbcut(iter)=no;

variables
lb  'the objective function of main problem (lower bound of original problem)'
;

parameter
vv(i,iter)
ww(i,j,iter);

equations
optimalitycut(iter)  'benders optimality cutting plane'
feasibilitycut(iter) 'benders feasibility cutting plane'
;

optimalitycut(optcut) .. lb=g=sum(j,y(j)*f(j))+sum(i,vv(i,optcut))-sum((i,j),ww(i,j,optcut)*y(j));
feasibilitycut(unbcut) .. sum(i,vv(i,unbcut))-sum((i,j),ww(i,j,unbcut)*y(j))=l=0;

model master /optimalitycut, feasibilitycut/;

*///////////////////////////// benders algorithm //////////////////////////
parameter results(iter,*);
parameter nonconverged /yes/;
scalar  unbounded  /1.0e9/;
scalar uperbound /+inf/;
scalar lowerbound /-inf/;

loop(iter$(nonconverged),

solve  dualsubproblem maximizing zdsp using lp;
if (zdsp.l<unbounded,
uperbound=zdsp.l;
optcut(iter)=yes;
else
solve  extremeray maximizing zdsp using lp;
unbcut(iter)=yes;
);

vv(i,iter)=v.l(i);
ww(i,j,iter)=w.l(i,j);

option optcr=0;
solve master minimizing lb using mip;
display lb.l;

abort$(master.modelstat<>1) "original model not solved optimality";
*abort$(dualsubproblem.modelstat=4)
lowerbound=lb.l;
results(iter,'lb')=lowerbound;
results(iter,'ub')=uperbound;
nonconverged$((uperbound-lowerbound)<0.1)=no;
);

display  uperbound, lowerbound;

