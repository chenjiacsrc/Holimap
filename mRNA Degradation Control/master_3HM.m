function f = master_3HM(t,x,para,input)

% -------------------------------------------------------------------------
% model parameters
h = para(1); v = para(2); alpha = para(3); d = para(4);
rhou = para(5); rhob = para(6);
sigmau = para(7); sigmab = para(8);
u = para(9); a = para(10); b = para(11);
N = para(12); M = para(13);

% -------------------------------------------------------------------------
% derived parameters
inputg0 = input(:,01); inputg1 = input(:,02);
inputm0 = input(:,03); inputm1 = input(:,04);
inputp0 = input(:,05); inputmn = input(:,06);
inputmp0 = input(:,07); time = input(:,08); 

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
temp = abs(time-t);
ind = find(temp==min(temp)); ind = min(ind);
g0 = inputg0(ind)+eps; g1 = inputg1(ind)+eps;
m0 = inputm0(ind)+eps; m1 = inputm1(ind)+eps;
p0 = inputp0(ind)+eps; mn = inputmn(ind)+eps;
mp0 = inputmp0(ind)+eps; m = m0+m1;
vbar = v+alpha*mn/m;
mat = zeros(2,2);
mat(1,:) = [g1,-g0];
mat(2,:) = [m1,-m0];
aim = zeros(2,1);
aim(1) = sigmau*g1-sigmab*p0;
aim(2) = sigmau*m1-sigmab*mp0;
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
    Q(i+1,i) = i*vbar;
    Q(i+1+M,i+M) = i*vbar;
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