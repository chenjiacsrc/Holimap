function f = master_2HM(t,x,para,input)

% -------------------------------------------------------------------------
% model parameters
h1 = para(1); h2 = para(2); h3 = para(3);
d1 = para(4); d2 = para(5); d3 = para(6);
rhou1 = para(7); rhob1 = para(8);
rhou2 = para(9); rhob2 = para(10);
rhou3 = para(11); rhob3 = para(12);
sigmau1 = para(13); sigmab1 = para(14);
sigmau2 = para(15); sigmab2 = para(16); 
sigmau3 = para(17); sigmab3 = para(18);
M = para(19); N = para(20); P = para(21);

% -------------------------------------------------------------------------
% derived parameters
inputg0 = input(:,01); inputg1 = input(:,02); inputg11 = input(:,03);
inputm0 = input(:,04); inputm1 = input(:,05);
inputm10 = input(:,06); inputmm10 = input(:,07); inputmmm10 = input(:,08);
inputp0 = input(:,09); inputpp0 = input(:,10); inputppp0 = input(:,11);
inputmp0 = input(:,12); inputmpp0 = input(:,13); inputmppp0 = input(:,14);
time = input(:,15);

% -------------------------------------------------------------------------
% find the optimal time
eps = 1e-5; numtime = 50;
temp = abs(time-t);
ind = find(temp==min(temp)); ind = min(ind);
g0 = inputg0(ind)+eps; g1 = inputg1(ind)+eps; g11 = inputg11(ind)+eps;
m0 = inputm0(ind)+eps; m1 = inputm1(ind)+eps;
m10 = inputm10(ind)+eps; mm10 = inputmm10(ind)+eps; mmm10 = inputmmm10(ind)+eps;
p0 = inputp0(ind)+eps; pp0 = inputpp0(ind)+eps; ppp0 = inputppp0(ind)+eps;
mp0 = inputmp0(ind)+eps; mpp0 = inputmpp0(ind)+eps; mppp0 = inputmppp0(ind)+eps;
mat = zeros(2,2);
mat(1,:) = [g1,-g0];
mat(2,:) = [m1,-m0];
aim = zeros(2,1);
if h1 == 1 && h3 == 1
    aim(1) = sigmau1*g1-sigmab1*p0;
    aim(2) = sigmau1*m1-sigmab1*mp0-h1*sigmau2*g11+h1*sigmab2*m10;
elseif h1 == 1 && h3 == 2
    aim(1) = sigmau1*g1-sigmab1*pp0;
    aim(2) = sigmau1*m1-sigmab1*mpp0-h1*sigmau2*g11+h1*sigmab2*m10;
elseif h1 == 1 && h3 == 3
    aim(1) = sigmau1*g1-sigmab1*ppp0;
    aim(2) = sigmau1*m1-sigmab1*mppp0-h1*sigmau2*g11+h1*sigmab2*m10;
elseif h1 == 2 && h3 == 1
    aim(1) = sigmau1*g1-sigmab1*p0;
    aim(2) = sigmau1*m1-sigmab1*mp0-h1*sigmau2*g11+h1*sigmab2*mm10;
elseif h1 == 2 && h3 == 2
    aim(1) = sigmau1*g1-sigmab1*pp0;
    aim(2) = sigmau1*m1-sigmab1*mpp0-h1*sigmau2*g11+h1*sigmab2*mm10;
elseif h1 == 2 && h3 == 3
    aim(1) = sigmau1*g1-sigmab1*ppp0;
    aim(2) = sigmau1*m1-sigmab1*mppp0-h1*sigmau2*g11+h1*sigmab2*mm10;
elseif h1 == 3 && h3 == 1
    aim(1) = sigmau1*g1-sigmab1*p0;
    aim(2) = sigmau1*m1-sigmab1*mp0-h1*sigmau2*g11+h1*sigmab2*mmm10;
elseif h1 == 3 && h3 == 2
    aim(1) = sigmau1*g1-sigmab1*pp0;
    aim(2) = sigmau1*m1-sigmab1*mpp0-h1*sigmau2*g11+h1*sigmab2*mmm10;
elseif h1 == 3 && h3 == 3
    aim(1) = sigmau1*g1-sigmab1*ppp0;
    aim(2) = sigmau1*m1-sigmab1*mppp0-h1*sigmau2*g11+h1*sigmab2*mmm10;
end
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
    Q(i+1,i) = i*d1;
    Q(i+1+M,i+M) = i*d1;
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