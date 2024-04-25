function f = master_2HM(t,x,para,input)

% -------------------------------------------------------------------------
% model parameters
h = para(1); d = para(2); rho1 = para(3); rho0 = para(4);
alpha1 = para(5); alpha0 = para(6); sigma = para(7); M = para(8);

% -------------------------------------------------------------------------
% derived parameters
inputg0 = input(:,01); inputg1 = input(:,02);
inputm0 = input(:,03); inputm1 = input(:,04); inputmm0 = input(:,05);
inputq0 = input(:,06); inputmq0 = input(:,07); time = input(:,08);

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
temp = abs(time-t);
ind = find(temp==min(temp)); ind = min(ind);
g0 = inputg0(ind)+eps; g1 = inputg1(ind)+eps;
m0 = inputm0(ind)+eps; m1 = inputm1(ind)+eps; mm0 = inputmm0(ind)+eps;
q0 = inputq0(ind)+eps; mq0 = inputmq0(ind)+eps;
mat = zeros(2,2);
mat(1,:) = [g1,-g0];
mat(2,:) = [m1,-m0];
aim = zeros(2,1);
aim(1) = alpha1*g1-alpha0*g0-sigma*(m0+q0);
aim(2) = alpha1*m1-alpha0*m0-sigma*(mm0+mq0);
temp = inv(mat)*aim;
sigma1bar = temp(1);
sigma0bar = temp(2);

% -------------------------------------------------------------------------
% generator matrix
num = 2*M;
Q = zeros(num);
one = ones(num,1);
for i = 1:M-1
    Q(i,i+1) = rho0;
    Q(i+M,i+1+M) = rho1;
    Q(i+1,i) = i*d;
    Q(i+1+M,i+M) = i*d;
end
for i = 1:M
    Q(i,i+M) = sigma0bar;
    Q(i+M,i) = sigma1bar;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;