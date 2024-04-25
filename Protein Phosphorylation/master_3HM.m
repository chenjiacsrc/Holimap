function f = master_3HM(t,x,para,input)

% -------------------------------------------------------------------------
% model parameters
h = para(1); d0 = para(2); d1 = para(3); d2 = para(4);
rhou = para(5); rhob = para(6);
sigmau = para(7); sigmab = para(8);
a1 = para(9); a2 = para(10); a3 = para(11); a4 = para(12);
b1 = para(13); b2 = para(14); b3 = para(15); b4 = para(16);
c1 = para(17); c2 = para(18); c3 = para(19); c4 = para(20);
N1 = para(21); N2 = para(22); N3 = para(23); N4 = para(24); M = para(25);

% -------------------------------------------------------------------------
% derived parameters
inputg0 = input(:,01); inputg1 = input(:,02);
inputm1 = input(:,03); inputm2 = input(:,04);
inputn0m1 = input(:,05); inputn00 = input(:,06); inputn01 = input(:,07); 
inputn20 = input(:,08); inputn0n20 = input(:,09); time = input(:,10);

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
temp = abs(time-t);
ind = find(temp==min(temp)); ind = min(ind);
g0 = inputg0(ind)+eps; g1 = inputg1(ind)+eps;
m1 = inputm1(ind)+eps; m2 = inputm2(ind)+eps;
n0m1 = inputn0m1(ind)+eps; n00 = inputn00(ind)+eps; n01 = inputn01(ind)+eps;
n20 = inputn20(ind)+eps; n0n20 = inputn0n20(ind)+eps;
n0 = n00+n01;
dbar = d0+a1*N1-a1*n0m1/n0-b1*m1/n0-c2*m2/n0;
mat = zeros(2,2);
mat(1,:) = [g1,-g0];
mat(2,:) = [n01,-n00];
aim = zeros(2,1);
aim(1) = sigmau*g1-sigmab*n20;
aim(2) = sigmau*n01-sigmab*n0n20;
temp = inv(mat)*aim;
sigmaubar = temp(1);
sigmabbar = temp(2);

% -------------------------------------------------------------------------
% generator matrix
num = 2*M;
Q = zeros(num);
one = ones(num,1);
for i = 1:M-1
    Q(i,i+1) = rhou;
    Q(i+M,i+1+M) = rhob;
    Q(i+1,i) = i*dbar;
    Q(i+1+M,i+M) = i*dbar;
end
for i = 1:M
    Q(i,i+M) = sigmabbar;
    Q(i+M,i) = sigmaubar;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;