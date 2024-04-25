function f = master_3HM(t,x,para,input)

% -------------------------------------------------------------------------
% model parameters
h1 = para(1); h2 = para(2);
d1 = para(3); d2 = para(4);
rhou1 = para(5); rhob1 = para(6);
rhou2 = para(7); rhob2 = para(8);
sigmau1 = para(9); sigmab1 = para(10);
sigmau2 = para(11); sigmab2 = para(12);
alpha = para(13); beta = para(14);
a1 = para(15); a2 = para(16);
b1 = para(17); b2 = para(18); M = para(19);

% -------------------------------------------------------------------------
% derived parameters
inputm = input(:,01); inputmr = input(:,02);
inputc1 = input(:,03); time = input(:,04);

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
temp = abs(time-t);
ind = find(temp==min(temp)); ind = min(ind);
m = inputm(ind)+eps; mr = inputmr(ind)+eps; c1 = inputc1(ind)+eps;
dbar = d1+alpha*mr/m-(beta+a2)*c1/m;

% -------------------------------------------------------------------------
% generator matrix
num = 2*M;
Q = zeros(num);
one = ones(num,1);
for i = 1:M-1
    Q(i,i+1) = rhou1;
    Q(i+M,i+1+M) = rhob1;
    Q(i+1,i) = i*dbar;
    Q(i+1+M,i+M) = i*dbar;
end
for i = 1:M
    Q(i,i+M) = sigmab1;
    Q(i+M,i) = sigmau1;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;