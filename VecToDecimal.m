function y = VecToDecimal(x,base)

z = repmat(base.^(size(x,2)-1:-1:0),size(x,1),1);
y = sum(x.*z,2);

end

