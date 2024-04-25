tic
warning off all

% -------------------------------------------------------------------------
% model parameters
h1 = 1; h2 = 1; d1 = 1; d2 = 1;
rhou1 = 0.5*d1; rhob1 = 3*d1; rhou2 = 0*d2; rhob2 = 3*d2;
sigmau1 = 1; sigmau2 = sigmau1; L1 = 30; L2 = 1;
sigmab1 = L1*sigmau1; sigmab2 = L2*sigmau2;
alpha = 0.1; beta = 1e-3; a1 = 2; a2 = 0.5; b1 = 2; b2 = 0.5;

% -------------------------------------------------------------------------
% derived parameters
K = 30; T = 10; Tc = T;
M = 121; maxm = M;
rhou1 = rhou1*K; rhob1 = rhob1*K;
rhou2 = rhou2*K; rhob2 = rhob2*K;
sigmab1 = sigmab1/K^h2; sigmab2 = sigmab2/K^h1; alpha = alpha/K;
para = [h1,h2,d1,d2,rhou1,rhob1,rhou2,rhob2,sigmau1,sigmab1,...
    sigmau2,sigmab2,alpha,beta,a1,a2,b1,b2,M];

% -------------------------------------------------------------------------
% simulate trajectories
step = 2e-2/K/sigmau1; num = 2e3;
time = 0:step:T; len = length(time);
x1 = zeros(num,len); x2 = zeros(num,len); x3 = zeros(num,len); x4 = zeros(num,len);
s1 = zeros(num,len); s2 = zeros(num,len);
for k = 1:num
    % initial value
    x1(k,1) = 0; x2(k,1) = 0; x3(k,1) = 0; x4(k,1) = 0;
    s1(k,1) = 1; s2(k,1) = 0;
    
    % evolution
    for i = 1:len-1
        rho1 = rhou1*(s1(k,i)==0)+rhob1*(s1(k,i)==1);
        rho2 = rhou2*(s2(k,i)==0)+rhob2*(s2(k,i)==1);
        rhotot = rho1+rho2;
        lambda1 = alpha*x1(k,i)*x2(k,i)+beta*x3(k,i);
        lambda2 = alpha*x2(k,i)*x3(k,i)+beta*x4(k,i);
        lambdatot = lambda1+lambda2;
        dtot = d1*x1(k,i)+d2*x2(k,i)+(a1+a2)*x3(k,i)+(b1+b2)*x4(k,i);
        sigma1 = sigmab1*(s1(k,i)==0)+sigmau1*(s1(k,i)==1);
        sigma2 = sigmab2*(s2(k,i)==0)+sigmau2*(s2(k,i)==1);
        sigmatot = sigma1+sigma2;
        tot = rhotot+lambdatot+dtot+sigmatot;
        if rand < step*tot
            temp = rand;
            if temp < rho1/tot
                x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < rhotot/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+1; x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+alpha*x1(k,i)*x2(k,i))/tot
                x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i)-1; x3(k,i+1) = x3(k,i)+1; x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambda1)/tot
                x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i)+1; x3(k,i+1) = x3(k,i)-1; x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambda1+alpha*x2(k,i)*x3(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-1; x3(k,i+1) = x3(k,i)-1; x4(k,i+1) = x4(k,i)+1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambdatot)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+1; x3(k,i+1) = x3(k,i)+1; x4(k,i+1) = x4(k,i)-1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambdatot+d1*x1(k,i))/tot
                x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambdatot+d1*x1(k,i)+d2*x2(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-1; x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambdatot+d1*x1(k,i)+d2*x2(k,i)+a1*x3(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+1; x3(k,i+1) = x3(k,i)-1; x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambdatot+d1*x1(k,i)+d2*x2(k,i)+(a1+a2)*x3(k,i))/tot
                x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i)-1; x4(k,i+1) = x4(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambdatot+d1*x1(k,i)+d2*x2(k,i)+(a1+a2)*x3(k,i)+b1*x4(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+2; x3(k,i+1) = x3(k,i)-1; x4(k,i+1) = x4(k,i)-1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);                   
            elseif temp < (rhotot+lambdatot+dtot)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i)+1; x4(k,i+1) = x4(k,i)-1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
            elseif temp < (rhotot+lambdatot+dtot+sigma1)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i); x4(k,i+1) = x4(k,i);
                s1(k,i+1) = 1-s1(k,i); s2(k,i+1) = s2(k,i);
            else
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = 1-s2(k,i);
            end
        else
            x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i);
            s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i);
        end
    end
end

% -------------------------------------------------------------------------
% downsample
numpoint = 200; itv = floor(T/step/numpoint);
time = time(1:itv:end)'; len = length(time);
x1 = x1(:,1:itv:end); x2 = x2(:,1:itv:end);
x3 = x3(:,1:itv:end); x4 = x4(:,1:itv:end);
s1 = s1(:,1:itv:end); s2 = s2(:,1:itv:end);

% -------------------------------------------------------------------------
% calculate moments
m = sum(x1)'/num;
mr = sum(x1.*x2)'/num;
c1 = sum(x3)'/num;
input = [m,mr,c1,time];

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
% ymax = 0.033; ylim([0,ymax]); set(gca,'ytick',[0,ymax/3,ymax*2/3,ymax]);
patch([hori fliplr(hori)], [vert 1e-4*ones(size(vert))],[0.9 0.9 0.9],'LineStyle','none'); hold on

% -------------------------------------------------------------------------
% predict distributions using Holimap
birthhm = zeros(1,2*M); birthhm(M+1) = 1;
[~,sol] = ode45(@master_3HM,[0,Tc],birthhm,[],para,input);
disthm = sol(end,1:M)+sol(end,M+1:2*M);
disthm = max(disthm,zeros(1,M));
plot(0:maxm-1,disthm,'b','Linewidth',2); hold on

toc