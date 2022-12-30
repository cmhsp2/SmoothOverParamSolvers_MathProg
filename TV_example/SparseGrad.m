function [T1,T2] = SparseGrad(n1,n2)

e = ones(n1,1);
e2 = ones(n2,1);

D1 =  spdiags([ -1*e e ],0:1,n1,n1);
D1(n1,n1)=0;
D2 =  spdiags([ -1*e2 e2 ],0:1,n2,n2);
D2(n2,n2)=0;

E1 =  speye(n1);
E2 =  speye(n2);

T1 = kron(D2,E1);
T2 = kron(E2,D1);
end