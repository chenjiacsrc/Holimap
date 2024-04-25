function f = master_LMA(para,input,birth,T)

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
tps = linspace(0.2*T,T,numtime);
sigmab1bar = 0;
for i = 1:numtime
    tp = tps(i);
    temp = abs(time-tp);
    ind = find(temp==min(temp)); ind = min(ind);
    g0 = inputg0(ind)+eps;
    p0 = inputp0(ind)+eps;
    pp0 = inputpp0(ind)+eps;
    ppp0 = inputppp0(ind)+eps;
    if h3 == 1
        sigmab1bar = sigmab1bar+sigmab1*p0/g0;
    elseif h3 == 2
        sigmab1bar = sigmab1bar+sigmab1*pp0/g0;
    elseif h3 == 3
        sigmab1bar = sigmab1bar+sigmab1*ppp0/g0;
    end
end
sigmab1bar = sigmab1bar/numtime;

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
    Q(i+M,i) = sigmau1;
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = birth*expm(T*Q);