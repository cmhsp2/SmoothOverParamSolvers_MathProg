
function tmp = groupthresh_TV(z,tau)

nm = sqrt(sum(z(1:end/2,:).^2 + z(end/2+1:end,:).^2,2)) ;
tmp = max( nm - tau,0);
INV = 1./nm;
INV(isinf(INV)) = 0;
tmp = repmat(tmp.*INV,2,1).*z;

end