clear
clc

   
M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];

Slist = [[0; 0;  1;  4; 0;    0], ...
       [0; 0;  0;  0; 1;    0], ...
       [0; 0; -1; -6; 0; -0.1]];
thetalist =[pi / 2; 3; pi];

Blist=SlistToBlist(Slist,M);
Jb_1 = JacobianBody(Blist, thetalist);


T = FKinSpace(M, Slist, thetalist);


Js = JacobianSpace(Slist, thetalist);

Jb_2=Adjoint(TransInv(T))*Js;