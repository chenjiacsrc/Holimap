tic
warning off all

% -------------------------------------------------------------------------
% model parameters
h = 1; d0 = 1; d1 = 0.1; d2 = 3;
rhou = 1; rhob = 5; sigmau = 1; sigmab = 5*sigmau;
a1 = 1; a2 = 0.1; a3 = a1; a4 = a2;
b1 = 0.1; b2 = 1; b3 = b1; b4 = b2;
c1 = 1; c2 = 0.1; c3 = c1; c4 = c2;
N1 = 100; N2 = N1; N3 = N1; N4 = N1;

% -------------------------------------------------------------------------
% derived parameters
K = 30; T = 10; Tc = T;
M = 121; maxm = M;
rhou = rhou*K; rhob = rhob*K;
sigmab = sigmab/K^h; a1 = a1/K; a2 = a2/K; a3 = a3/K; a4 = a4/K;
para = [h,d0,d1,d2,rhou,rhob,sigmau,sigmab,a1,a2,a3,a4,b1,b2,b3,b4,...
    c1,c2,c3,c4,N1,N2,N3,N4,M];

% -------------------------------------------------------------------------
% simulate trajectories
step = 2e-2/K/sigmau; num = 2e3;
time = 0:step:T; len = length(time);
x0 = zeros(num,len); x1 = zeros(num,len); x2 = zeros(num,len); s = zeros(num,len);
y1 = zeros(num,len); y2 = zeros(num,len); y3 = zeros(num,len); y4 = zeros(num,len);
for k = 1:num
    % initial value
    x0(k,1) = 0; x1(k,1) = 0; x2(k,1) = 0; s(k,1) = 1;
    y1(k,1) = 0; y2(k,1) = 0; y3(k,1) = 0; y4(k,1) = 0;
    
    % evolution
    for i = 1:len-1
        rhotot = rhou*(s(k,i)==0)+rhob*(s(k,i)==1);
        dtot = d0*x0(k,i)+d1*x1(k,i)+d2*x2(k,i);
        atot = a1*x0(k,i)*(N1-y1(k,i))+a2*x1(k,i)*(N2-y2(k,i))+a3*x1(k,i)*(N3-y3(k,i))+a4*x2(k,i)*(N4-y4(k,i));
        btot = b1*y1(k,i)+b2*y2(k,i)+b3*y3(k,i)+b4*y4(k,i);
        ctot = c1*y1(k,i)+c2*y2(k,i)+c3*y3(k,i)+c4*y4(k,i);
        prod = 1;
        for ii = 1:h
            prod = prod*(x2(k,i)-ii+1);
        end
        sigmatot = sigmab*prod*(s(k,i)==0)+sigmau*(s(k,i)==1);
        tot = rhotot+dtot+atot+btot+ctot+sigmatot;
        if rand < step*tot
            temp = rand;
            if temp < rhotot/tot
                x0(k,i+1) = x0(k,i)+1; x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+d0*x0(k,i))/tot
                x0(k,i+1) = x0(k,i)-1; x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+d0*x0(k,i)+d1*x1(k,i))/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot)/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-1; s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+a1*x0(k,i)*(N1-y1(k,i)))/tot
                x0(k,i+1) = x0(k,i)-1; x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i)+1; y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+a1*x0(k,i)*(N1-y1(k,i))+a2*x1(k,i)*(N2-y2(k,i)))/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i)+1; y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+a1*x0(k,i)*(N1-y1(k,i))+a2*x1(k,i)*(N2-y2(k,i))+a3*x1(k,i)*(N3-y3(k,i)))/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i)+1; y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+atot)/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-1; s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i)+1;
            elseif temp < (rhotot+dtot+atot+b1*y1(k,i))/tot
                x0(k,i+1) = x0(k,i)+1; x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i)-1; y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+atot+b1*y1(k,i)+b2*y2(k,i))/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i)-1; y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+atot+b1*y1(k,i)+b2*y2(k,i)+b3*y3(k,i))/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i)-1; y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+atot+btot)/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+1; s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i)-1;
            elseif temp < (rhotot+dtot+atot+btot+c1*y1(k,i))/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i)-1; y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+atot+btot+c1*y1(k,i)+c2*y2(k,i))/tot
                x0(k,i+1) = x0(k,i)+1; x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i)-1; y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+atot+btot+c1*y1(k,i)+c2*y2(k,i)+c3*y3(k,i))/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+1; s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i)-1; y4(k,i+1) = y4(k,i);
            elseif temp < (rhotot+dtot+atot+btot+ctot)/tot
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i)-1;
            else
                x0(k,i+1) = x0(k,i)-h*(s(k,i)==0)+h*(s(k,i)==1); x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); s(k,i+1) = 1-s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
            end
        else
                x0(k,i+1) = x0(k,i); x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); s(k,i+1) = s(k,i);
                y1(k,i+1) = y1(k,i); y2(k,i+1) = y2(k,i); y3(k,i+1) = y3(k,i); y4(k,i+1) = y4(k,i);
        end
    end
end

% -------------------------------------------------------------------------
% downsample
numpoint = 200; itv = floor(T/step/numpoint);
time = time(1:itv:end)'; len = length(time);
x0 = x0(:,1:itv:end); x1 = x1(:,1:itv:end); x2 = x2(:,1:itv:end); s = s(:,1:itv:end);
y1 = y1(:,1:itv:end); y2 = y2(:,1:itv:end); y3 = y3(:,1:itv:end); y4 = y4(:,1:itv:end);

% -------------------------------------------------------------------------
% calculate moments
g0 = sum((s(:,:)==0))'/num;
g1 = sum((s(:,:)==1))'/num;
m1 = sum(y1)'/num;
m2 = sum(y2)'/num;
n0m1 = sum(x0.*y1)'/num;
n00 = sum(x0.*(s(:,:)==0))'/num;
n01 = sum(x0.*(s(:,:)==1))'/num;
n20 = sum(x2.*(s(:,:)==0))'/num;
n0n20 = sum(x0.*x2.*(s(:,:)==0))'/num;
input = [g0,g1,m1,m2,n0m1,n00,n01,n20,n0n20,time];

% -------------------------------------------------------------------------
% simulate distributions using SSA
backnum = 0;
temp = abs(time-Tc);
ind = find(temp==min(temp)); ind = min(ind);
distm = zeros(1,M);
for i = 1:M
    distm(i) = sum(sum(x0(:,ind-backnum:ind+backnum)==i-1))/num/(2*backnum+1);
end
hori = 0:maxm-1; vert = distm(1:maxm);
figure; plot(hori,vert,'color',0.6*ones(3,1)); box off; hold on
xmax = 120; xlim([0,xmax]); set(gca,'xtick',[0,xmax/3,xmax*2/3,xmax]);
ymax = 0.03; ylim([0,ymax]); set(gca,'ytick',[0,ymax/3,ymax*2/3,ymax]);
% ymax = 0.048; ylim([0,ymax]); set(gca,'ytick',[0,ymax/3,ymax*2/3,ymax]);
patch([hori fliplr(hori)], [vert 1e-4*ones(size(vert))],[0.9 0.9 0.9],'LineStyle','none'); hold on

% -------------------------------------------------------------------------
% predict distributions using Holimap
birthhm = zeros(1,2*M); birthhm(M+1) = 1;
[~,sol] = ode45(@master_3HM,[0,Tc],birthhm,[],para,input);
disthm = sol(end,1:M)+sol(end,M+1:2*M);
disthm = max(disthm,zeros(1,M));
plot(0:maxm-1,disthm,'b','Linewidth',2); hold on

toc