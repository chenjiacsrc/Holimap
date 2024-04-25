tic
warning off all

% -------------------------------------------------------------------------
% model parameters
h = 1; d = 1; rho1 = 81; rho0 = 5.4;
alpha1 = 1; alpha0 = 0.5; sigma = 0.02;

% -------------------------------------------------------------------------
% derived parameters
T = 10; Tc = T;
M = 121; maxm = M;
para = [h,d,rho1,rho0,alpha1,alpha0,sigma,M];

% -------------------------------------------------------------------------
% simulate trajectories
step = 2e-3; num = 2e3;
time = 0:step:T; len = length(time);
x1 = zeros(num,len); x2 = zeros(num,len); x3 = zeros(num,len); x4 = zeros(num,len);
x5 = zeros(num,len); x6 = zeros(num,len); x7 = zeros(num,len); x8 = zeros(num,len);
s1 = zeros(num,len); s2 = zeros(num,len); s3 = zeros(num,len); s4 = zeros(num,len);
s5 = zeros(num,len); s6 = zeros(num,len); s7 = zeros(num,len); s8 = zeros(num,len);
for k = 1:num
    % initial value
    x1(k,1) = 0; x2(k,1) = 0; x3(k,1) = 0; x4(k,1) = 0;
    x5(k,1) = 0; x6(k,1) = 0; x7(k,1) = 0; x8(k,1) = 0;
    s1(k,1) = 1; s2(k,1) = 0; s3(k,1) = 0; s4(k,1) = 0;
    s5(k,1) = 0; s6(k,1) = 0; s7(k,1) = 0; s8(k,1) = 0;
    
    % evolution
    for i = 1:len-1
        r1 = rho0*(s1(k,i)==0)+rho1*(s1(k,i)==1);
        r2 = rho0*(s2(k,i)==0)+rho1*(s2(k,i)==1);
        r3 = rho0*(s3(k,i)==0)+rho1*(s3(k,i)==1);
        r4 = rho0*(s4(k,i)==0)+rho1*(s4(k,i)==1);
        r5 = rho0*(s5(k,i)==0)+rho1*(s5(k,i)==1);
        r6 = rho0*(s6(k,i)==0)+rho1*(s6(k,i)==1);
        r7 = rho0*(s7(k,i)==0)+rho1*(s7(k,i)==1);
        r8 = rho0*(s8(k,i)==0)+rho1*(s8(k,i)==1);
        rhotot = r1+r2+r3+r4+r5+r6+r7+r8;
        dtot = d*(x1(k,i)+x2(k,i)+x3(k,i)+x4(k,i)+x5(k,i)+x6(k,i)+x7(k,i)+x8(k,i));
        sigma1 = (alpha0+sigma*(x1(k,i)+x6(k,i)))*(s1(k,i)==0)+(alpha1)*(s1(k,i)==1);
        sigma2 = (alpha0)*(s2(k,i)==0)+(alpha1+sigma*x1(k,i))*(s2(k,i)==1);
        sigma3 = (alpha0+sigma*x4(k,i))*(s3(k,i)==0)+(alpha1+sigma*x2(k,i))*(s3(k,i)==1);
        sigma4 = (alpha0)*(s4(k,i)==0)+(alpha1)*(s4(k,i)==1);
        sigma5 = (alpha0+sigma*x3(k,i))*(s5(k,i)==0)+(alpha1+sigma*(x1(k,i)+x5(k,i)))*(s5(k,i)==1);
        sigma6 = (alpha0+sigma*x5(k,i))*(s6(k,i)==0)+(alpha1)*(s6(k,i)==1);
        sigma7 = (alpha0+sigma*x3(k,i))*(s7(k,i)==0)+(alpha1+sigma*x6(k,i))*(s7(k,i)==1);
        sigma8 = (alpha0+sigma*x7(k,i))*(s8(k,i)==0)+(alpha1)*(s8(k,i)==1);
        sigmatot = sigma1+sigma2+sigma3+sigma4+sigma5+sigma6+sigma7+sigma8;
        tot = rhotot+dtot+sigmatot;
        if rand < step*tot
            temp = rand;
            if temp < r1/tot
                x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (r1+r2)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+1; x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (r1+r2+r3)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i)+1; x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (r1+r2+r3+r4)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i)+1;
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (r1+r2+r3+r4+r5)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i)+1; x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (r1+r2+r3+r4+r5+r6)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i)+1; x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (r1+r2+r3+r4+r5+r6+r7)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i)+1; x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < rhotot/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i)+1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+d*x1(k,i))/tot
                x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+d*x1(k,i)+d*x2(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-1; x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+d*x1(k,i)+d*x2(k,i)+d*x3(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i)-1; x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+d*x1(k,i)+d*x2(k,i)+d*x3(k,i)+d*x4(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i)-1;
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+d*x1(k,i)+d*x2(k,i)+d*x3(k,i)+d*x4(k,i)+d*x5(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i)-1; x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+d*x1(k,i)+d*x2(k,i)+d*x3(k,i)+d*x4(k,i)+d*x5(k,i)+d*x6(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i)-1; x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+d*x1(k,i)+d*x2(k,i)+d*x3(k,i)+d*x4(k,i)+d*x5(k,i)+d*x6(k,i)+d*x7(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i)-1; x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i)-1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot+sigma1)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = 1-s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot+sigma1+sigma2)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = 1-s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot+sigma1+sigma2+sigma3)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = 1-s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot+sigma1+sigma2+sigma3+sigma4)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = 1-s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot+sigma1+sigma2+sigma3+sigma4+sigma5)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = 1-s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot+sigma1+sigma2+sigma3+sigma4+sigma5+sigma6)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = 1-s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = s8(k,i);
            elseif temp < (rhotot+dtot+sigma1+sigma2+sigma3+sigma4+sigma5+sigma6+sigma7)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = 1-s7(k,i); s8(k,i+1) = s8(k,i);
            else
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                x5(k,i+1) = x5(k,i); x6(k,i+1) = x6(k,i); x7(k,i+1) = x7(k,i); x8(k,i+1) = x8(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
                s5(k,i+1) = s5(k,i); s6(k,i+1) = s6(k,i); s7(k,i+1) = s7(k,i); s8(k,i+1) = 1-s8(k,i);
            end
        else
            x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
            s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i); s4(k,i+1) = s4(k,i);
        end
    end
end

% -------------------------------------------------------------------------
% downsample
numpoint = 200; itv = floor(T/step/numpoint);
time = time(1:itv:end)'; len = length(time);
x1 = x1(:,1:itv:end); x2 = x2(:,1:itv:end); x3 = x3(:,1:itv:end); x4 = x4(:,1:itv:end);
x5 = x5(:,1:itv:end); x6 = x6(:,1:itv:end); x7 = x7(:,1:itv:end); x8 = x8(:,1:itv:end);
s1 = s1(:,1:itv:end); s2 = s2(:,1:itv:end); s3 = s3(:,1:itv:end); s4 = s4(:,1:itv:end);
s5 = s5(:,1:itv:end); s6 = s6(:,1:itv:end); s7 = s7(:,1:itv:end); s8 = s8(:,1:itv:end);

% -------------------------------------------------------------------------
% calculate moments
g0 = sum((s1(:,:)==0))'/num;
g1 = sum((s1(:,:)==1))'/num;
g11 = sum((s1(:,:)==1).*(s2(:,:)==1))'/num;
m0 = sum(x1.*(s1(:,:)==0))'/num;
m1 = sum(x1.*(s1(:,:)==1))'/num;
mm0 = sum(x1.*(x1-1).*(s1(:,:)==0))'/num;
q0 = sum(x6.*(s1(:,:)==0))'/num;
mq0 = sum(x1.*x6.*(s1(:,:)==0))'/num;
input = [g0,g1,m0,m1,mm0,q0,mq0,time];

% -------------------------------------------------------------------------
% simulate distributions using SSA
backnum = 0;
temp = abs(time-Tc);
ind = find(temp==min(temp)); ind = min(ind);
distm = zeros(1,M);
for i = 1:M
    distm(i) = sum(sum(x1(:,ind-backnum:ind+backnum)==i-1))/num/(2*backnum+1);
end
hori = 0:maxm-1; vert = distm(1:maxm);
figure; plot(hori,vert,'color',0.6*ones(3,1)); box off; hold on
xmax = 120; xlim([0,xmax]); set(gca,'xtick',[0,xmax/3,xmax*2/3,xmax]);
ymax = 0.021; ylim([0,ymax]); set(gca,'ytick',[0,ymax/3,ymax*2/3,ymax]);
patch([hori fliplr(hori)], [vert 1.5e-4*ones(size(vert))],[0.9 0.9 0.9],'LineStyle','none'); hold on

% -------------------------------------------------------------------------
% predict distributions using Holimap
birthlma = zeros(1,2*M); birthlma(M+1) = 1;
[~,solhm] = ode45(@master_2HM,[0,Tc],birthlma,[],para,input);
disthm = solhm(end,1:M)+solhm(end,M+1:2*M);
disthm = max(disthm,zeros(1,M));
plot(0:maxm-1,disthm(1:maxm),'b','Linewidth',2); hold on

toc