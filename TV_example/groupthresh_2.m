
function tmp = groupthresh_2(z,tau)

nm = sqrt(sum(z.^2,2)) ;
tmp = max( nm - tau,0);
INV = 1./nm;
INV(isinf(INV)) = 0;
tmp = tmp.*INV.*z;

end