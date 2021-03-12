function [adj] = adjoint(A)
    a=A(1:3,4);ph=[0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
    adj=[A(1:3,1:3),ph*A(1:3,1:3);zeros(3),A(1:3,1:3)];
end