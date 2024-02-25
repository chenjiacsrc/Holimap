function f = master_LMA(t,x,para,input,momor,T)

% -------------------------------------------------------------------------
% model parameters
h = para(1); d = para(2); B = para(3); rhou = para(4); rhob = para(5);
sigmau = para(6); sigmab = para(7); N = para(8);

% -------------------------------------------------------------------------
% derived parameters
p = B/(B+1); q = 1-p; temp = momor;
time = input(1,:); inputg0 = input(2,:);
inputp0 = input(3,:); inputp1 = input(4,:);
inputpp0 = input(5,:); inputpp1 = input(6,:);

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
tps = linspace(0.2*T,T,numtime);
sigmabbar = 0;
for i = 1:numtime
    tp = tps(i);
    temp = abs(time-tp);
    ind = find(temp==min(temp)); ind = min(ind);
    g0 = inputg0(ind)+eps;
    p0 = inputp0(ind)+eps;
    pp0 = inputpp0(ind)+eps;
    if h == 1
        sigmabbar = sigmabbar+sigmab*p0/g0;
    else
        sigmabbar = sigmabbar+sigmab*pp0/g0;
    end
end
sigmabbar = sigmabbar/numtime;

% -------------------------------------------------------------------------
% generator matrix
num = 2*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+1+N,i+N) = i*d;
    for j = (i+1):N
        Q(i,j) = rhou*p^(j-i)*q;
        Q(i+N,j+N) = rhob*p^(j-i)*q;
    end
end
for i = 1:N
    Q(i,i+N) = sigmabbar;
    Q(i+N,i) = sigmau;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;