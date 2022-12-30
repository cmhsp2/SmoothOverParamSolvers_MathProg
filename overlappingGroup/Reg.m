%given a vector x and I a set of indices, return sum_{S in I} sqrt(|S|)* |x_S|_2

function R = Reg(x,I)
R = 0;
for i=1:length(I)
    R = R+ sqrt(length(I{i}))* norm(x(I{i}));
end

end
