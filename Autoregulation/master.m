function f = master(para)

% -------------------------------------------------------------------------
% model parameters
h = para(1); d = para(2); B = para(3); rhou = para(4); rhob = para(5);
sigmau = para(6); sigmab = para(7); N = para(8);

% -------------------------------------------------------------------------
% derived parameters
p = B/(B+1); q = 1-p;

% -------------------------------------------------------------------------
% generator matrix
num = 2*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    for j = (i+1):N
        Q(i,j) = rhou*p^(j-i)*q;
        Q(i+N,j+N) = rhob*p^(j-i)*q;
    end
    Q(i+1,i) = i*d;
    Q(i+1+N,i+N) = i*d;
end
for i = 1+h:N
    prop = prod(i-h:i-1);
    Q(i,i-h+N) = sigmab*prop;
    Q(i-h+N,i) = sigmau;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q;