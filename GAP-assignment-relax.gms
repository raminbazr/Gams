sets
i        task /1*3/
j        machine /1*2/
iter /iter1*iter100/;


table p(i,j)     profit
         1       2
1       10       8
2       13       14
3       9        10;

parameter d(j)   capacity of machine
/1 23,2 27/;

table w(i,j) amount of capacity of machine j used by task i
         1       2
1        10      9
2        14      11
3        18      13 ;

parameter u1(i) lagrangian multiplier;

u1('1')=0;
u1('2')=0;
u1('3')=0;

parameter results1(iter,*);
parameter results2(iter,*);
parameter results3(iter,*);

variable
z;
binary variable y(i,j);

equation
obj
cons2(j)

;

obj           .. z=e= sum((i,j),(p(i,j)-u1('1')-u1('2')-u1('3'))*y(i,j))+u1('1')+u1('2')+u1('3');
cons2(j)      .. sum(i,w(i,j)*y(i,j))=l= d(j);


model GAP /all/;

loop (iter$ (u1('1')<=1 and u1('2')<=1 and u1('3')<=1 ),

solve GAP using mip max z;
option mip=cplex;

results1(iter,'u1(1)')= u1('1');
results2(iter,'u1(2)')= u1('2');
results3(iter,'u1(3)')= u1('3');
results1(iter,'z')= z.l;
results2(iter,'z')= z.l;
results3(iter,'z')= z.l;

u1('1')= u1('1')+0.02;
u1('2')= u1('2')+0.02;
u1('3')= u1('3')+0.02;
);

display results1,results2,results3;

execute_unload 'Gap1.gdx'