function [gap] = gap_mat(A)
a=A(1:3,4);ph=[0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
gap=[eye(3),-ph;zeros(3),eye(3)];
end