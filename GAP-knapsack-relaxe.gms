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

parameter u2(j) lagrangian multiplier;

u2('1')=0;
u2('2')=0;

parameter results1(iter,*);
parameter results2(iter,*);

variable
z;
binary variable y(i,j);

equation
obj
cons1(i)

;

obj           .. z=e= sum((i,j),(p(i,j)-u2(j)*w(i,j))*y(i,j))+(u2('1')*23 + u2('2')*27);
cons1(i)      .. sum(j,y(i,j))=e=1;


model GAP /all/;

loop (iter$ (u2('1')<=1 and u2('2')<=1),

solve GAP using mip max z;
option mip=cplex;

results1(iter,'u2(1)')= u2('1');
results2(iter,'u2(2)')= u2('2');
results1(iter,'z')= z.l;
results2(iter,'z')= z.l;

u2('1')= u2('1')+0.02;
u2('2')= u2('2')+0.02;
);

display results1,results2;

execute_unload 'Gap.gdx'



