function f = master_LMA(t,x,para,input,T)

% -------------------------------------------------------------------------
% model parameters
h = para(end,1); d = para(end,2); B = para(end,3);
rhou = para(end,4); rhob = para(end,5);
sigmau = para(end,6); sigmab = para(end,7); N = para(end,8);

% -------------------------------------------------------------------------
% derived parameters
p = B/(B+1); q = 1-p;
temp = input; [len,num] = size(para); len = len-1;
input = para(1:len,1:num); time = para(1:len,end)';
momg0 = zeros(1,len);
momp0 = zeros(1,len); momp1 = zeros(1,len);
mompp0 = zeros(1,len); mompp1 = zeros(1,len);
momppp0 = zeros(1,len); momppp1 = zeros(1,len);
for k = 1:len
    momg0(k) = sum(input(k,1:N));
    for i = 0:N-1
        momp0(k) = momp0(k)+i*input(k,i+1);
        momp1(k) = momp1(k)+i*input(k,i+1+N);
        mompp0(k) = mompp0(k)+i*(i-1)*input(k,i+1);
        mompp1(k) = mompp1(k)+i*(i-1)*input(k,i+1+N);
        momppp0(k) = momppp0(k)+i*(i-1)*(i-2)*input(k,i+1);
        momppp1(k) = momppp1(k)+i*(i-1)*(i-2)*input(k,i+1+N);
    end
end
time = temp(1,:); inputg0 = temp(2,:);
inputp0 = temp(3,:); inputpp0 = temp(5,:);

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