function f = master_3HM(t,x,para,input)

% -------------------------------------------------------------------------
% model parameters
h1 = para(1); h2 = para(2);
d1 = para(3); d2 = para(4); alpha = para(5);
rhou1 = para(6); rhob1 = para(7);
rhou2 = para(8); rhob2 = para(9);
sigmau1 = para(10); sigmab1 = para(11);
sigmau2 = para(12); sigmab2 = para(13); M = para(14);

% -------------------------------------------------------------------------
% derived parameters
inputg0 = input(:,01); inputg1 = input(:,02);
inputm0 = input(:,03); inputm1 = input(:,04);
inputmm0 = input(:,05); inputmn1 = input(:,06); time = input(:,07);

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
temp = abs(time-t);
ind = find(temp==min(temp)); ind = min(ind);
g0 = inputg0(ind)+eps; g1 = inputg1(ind)+eps;
m0 = inputm0(ind)+eps; m1 = inputm1(ind)+eps;
mm0 = inputmm0(ind)+eps; mn1 = inputmn1(ind)+eps;
dbar = d1+alpha*mn1/m1;
mat = zeros(2,2);
mat(1,:) = [g1,-g0];
mat(2,:) = [m1,-m0];
aim = zeros(2,1);
aim(1) = sigmau1*g1-sigmab1*m0;
aim(2) = sigmau1*m1-sigmab1*mm0;
temp = inv(mat)*aim;
sigmau1bar = temp(1);
sigmab1bar = temp(2);

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
    Q(i,i+M) = sigmab1bar;
    Q(i+M,i) = sigmau1bar;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;