function f = master_4HM(t,x,para,input,momor)

% -------------------------------------------------------------------------
% model parameters
h = para(1); d = para(2); B = para(3); rhou = para(4); rhob = para(5);
sigmau = para(6); sigmab = para(7); N = para(8);

% -------------------------------------------------------------------------
% derived parameters
p = B/(B+1); q = 1-p; temp = input; input = momor;
time = input(1,:); inputg0 = input(2,:);
inputp0 = input(3,:); inputp1 = input(4,:);
inputpp0 = input(5,:); inputpp1 = input(6,:);
inputppp0 = input(7,:); inputppp1 = input(8,:);

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
temp = abs(time-t);
ind = find(temp==min(temp)); ind = min(ind);
g0 = inputg0(ind)+eps; g1 = 1-g0+eps;
p0 = inputp0(ind)+eps; p1 = inputp1(ind)+eps;
pp0 = inputpp0(ind)+eps; pp1 = inputpp1(ind)+eps;
ppp0 = inputppp0(ind)+eps; ppp1 = inputppp1(ind)+eps;

% ---------------------------------------------------------------------
% calculate the effective synthesis rates
mat = zeros(2,2);
mat(1,:) = [B*g0,B*g1];
mat(2,:) = [B*(p0+B*g0),B*(p1+B*g1)];
aim = zeros(2,1);
if h == 1
    aim(1) = rhou*B*g0+rhob*B*g1+sigmau*g1-sigmab*p0;
    aim(2) = rhou*B*(p0+B*g0)+rhob*B*(p1+B*g1)+sigmau*p1-sigmab*pp0;
else
    aim(1) = rhou*B*g0+rhob*B*g1+2*sigmau*g1-2*sigmab*pp0;
    aim(2) = rhou*B*(p0+B*g0)+rhob*B*(p1+B*g1)+sigmau*(2*p1+g1)-sigmab*(2*ppp0+pp0);
end
temp = inv(mat)*aim;
rhoubar = temp(1);
rhobbar = temp(2);

% ---------------------------------------------------------------------
% calculate the effective binding and unbinding rates
mat = zeros(2,2);
mat(1,:) = [g1,-g0];
mat(2,:) = [p1,-p0];
aim = zeros(2,1);
if h == 1
    aim(1) = sigmau*g1-sigmab*p0;
    aim(2) = sigmau*p1-sigmab*pp0+(rhobbar-rhob)*B*g1;
else
    aim(1) = sigmau*g1-sigmab*pp0;
    aim(2) = sigmau*p1-sigmab*ppp0+(rhobbar-rhob)*B*g1;
end
temp = inv(mat)*aim;
sigmaubar = temp(1);
sigmabbar = temp(2);

% -------------------------------------------------------------------------
% generator matrix
num = 2*N;
Q = zeros(num);
one = ones(num,1);
for i = 1:N-1
    Q(i+1,i) = i*d;
    Q(i+1+N,i+N) = i*d;
    for j = (i+1):N
        Q(i,j) = rhoubar*p^(j-i)*q;
        Q(i+N,j+N) = rhobbar*p^(j-i)*q;
    end
end
for i = 1:N
    Q(i,i+N) = sigmabbar;
    Q(i+N,i) = sigmaubar;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q'*x;