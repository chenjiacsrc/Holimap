tic
warning off all

% -------------------------------------------------------------------------
% model parameters
h = 1; v = 1; d = 0.2; alpha = 0.2;
rhou = 3.5; rhob = 0.7; sigmau = 2; sigmab = 30*sigmau;
u = 3; a = 1; b = 1; N = 100;

% -------------------------------------------------------------------------
% derived parameters
K = 30; T = 10; Tc = T;
M = 121; maxm = M;
rhou = rhou*K; rhob = rhob*K;
sigmab = sigmab/K^h; alpha = alpha/K;
para = [h,v,alpha,d,rhou,rhob,sigmau,sigmab,u,a,b,N,M];

% -------------------------------------------------------------------------
% simulate trajectories
step = 2e-2/K/sigmau; num = 2e3;
time = 0:step:T; len = length(time);
x1 = zeros(num,len); x2 = zeros(num,len);
x3 = zeros(num,len); s = zeros(num,len);
for k = 1:num
    % initial value
    x1(k,1) = 0; x2(k,1) = 0;
    x3(k,1) = 0; s(k,1) = 1;
    
    % evolution
    for i = 1:len-1
        rhotot = rhou*(s(k,i)==0)+rhob*(s(k,i)==1);
        dtot = v*x1(k,i)+d*x2(k,i)+alpha*x1(k,i)*x3(k,i);
        stot = a*(N-x3(k,i))+b*x3(k,i);
        prod = 1;
        for ii = 1:h
            prod = prod*(x2(k,i)-ii+1);
        end
        sigmatot = sigmab*prod*(s(k,i)==0)+sigmau*(s(k,i)==1);
        tot = rhotot+dtot+stot+sigmatot;
        if rand < step*tot
            temp = rand;
            if temp < rhotot/tot
                x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i);
                x3(k,i+1) = x3(k,i); s(k,i+1) = s(k,i);
            elseif temp < (rhotot+v*x1(k,i))/tot
                x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i);
                x3(k,i+1) = x3(k,i); s(k,i+1) = s(k,i);
            elseif temp < (rhotot+v*x1(k,i)+d*x2(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-1;
                x3(k,i+1) = x3(k,i); s(k,i+1) = s(k,i);
            elseif temp < (rhotot+dtot)/tot
                x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i);
                x3(k,i+1) = x3(k,i); s(k,i+1) = s(k,i);
            elseif temp < (rhotot+dtot+a*(N-x3(k,i)))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i);
                x3(k,i+1) = x3(k,i)+1; s(k,i+1) = s(k,i);
            elseif temp < (rhotot+dtot+stot)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i);
                x3(k,i+1) = x3(k,i)-1; s(k,i+1) = s(k,i);
            else
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-h*(s(k,i)==0)+h*(s(k,i)==1);
                x3(k,i+1) = x3(k,i); s(k,i+1) = 1-s(k,i);
            end
        else
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i);
                x3(k,i+1) = x3(k,i); s(k,i+1) = s(k,i);
        end
    end
end

% -------------------------------------------------------------------------
% downsample
numpoint = 200; itv = floor(T/step/numpoint);
time = time(1:itv:end)'; len = length(time);
x1 = x1(:,1:itv:end); x2 = x2(:,1:itv:end);
x3 = x3(:,1:itv:end); s = s(:,1:itv:end);

% -------------------------------------------------------------------------
% calculate moments
g0 = sum((s(:,:)==0))'/num;
g1 = sum((s(:,:)==1))'/num;
m0 = sum(x1.*(s(:,:)==0))'/num;
m1 = sum(x1.*(s(:,:)==1))'/num;
p0 = sum(x2.*(s(:,:)==0))'/num;
mn = sum(x1.*x3)'/num;
mp0 = sum(x1.*x2.*(s(:,:)==0))'/num;
input = [g0,g1,m0,m1,p0,mn,mp0,time];

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
ymax = 0.033; ylim([0,ymax]); set(gca,'ytick',[0,ymax/3,ymax*2/3,ymax]);
% ymax = 0.03; ylim([0,ymax]); set(gca,'ytick',[0,ymax/3,ymax*2/3,ymax]);
patch([hori fliplr(hori)], [vert 1e-4*ones(size(vert))],[0.9 0.9 0.9],'LineStyle','none'); hold on

% -------------------------------------------------------------------------
% predict distributions using Holimap
birthhm = zeros(1,2*M); birthhm(M+1) = 1;
[~,sol] = ode45(@master_3HM,[0,Tc],birthhm,[],para,input);
disthm = sol(end,1:M)+sol(end,M+1:2*M);
disthm = max(disthm,zeros(1,M));
plot(0:maxm-1,disthm,'b','Linewidth',2); hold on

toc