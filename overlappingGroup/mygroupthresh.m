function z2 = mygroupthresh(z,tau,idx)
z2 = zeros(size(z));
for i=1:length(idx)-1
    mid = idx(i):idx(i+1)-1;
    nm = norm(z(mid)) ;
    tmp = max( nm - tau,0);
    INV = 1./nm;
    INV(isinf(INV)) = 0;
    z2(mid) = (tmp.*INV).*z(mid);
    
end

end

