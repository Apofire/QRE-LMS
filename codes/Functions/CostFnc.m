function J = CostFnc(wts)

ITER = size(wts,2);
M    = size(wts,1);

J = zeros(ITER-1,1);
D = zeros(M,M);
W = zeros(M,1);

for i = 1:ITER-1
    for j = 1:M
        D(j,j) = 1/wts(j,i); 
        W(j)   = wts(j,i+1) - wts(j,i);
    end
    F = D\W;
    J(i) = norm(F);
end

end