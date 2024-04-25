tic
warning off all

% -------------------------------------------------------------------------
% model parameters
h1 = 3; h2 = 3; h3 = 3; d1 = 1; d2 = 1; d3 = 1;
rhou1 = 3*d1; rhob1 = 0.2*d1; rhou2 = 3*d2; rhob2 = 0*d2; rhou3 = 3*d3; rhob3 = 0*d3;
sigmau1 = 1.1; sigmau2 = sigmau1; sigmau3 = sigmau1; L1 = 10; L2 = L1; L3 = 1e3;
sigmab1 = L1*sigmau1; sigmab2 = L2*sigmau2; sigmab3 = L3*sigmau3;

% -------------------------------------------------------------------------
% derived parameters
K = 27; F = floor(K*0.5); T = 7; Tc = T;
M = 10*F; N = 10*F; P = 10*F; maxm = M; maxn = N; maxp = P;
rhou1 = rhou1*K; rhob1 = rhob1*K;
rhou2 = rhou2*K; rhob2 = rhob2*K;
rhou3 = rhou3*K; rhob3 = rhob3*K;
sigmab1 = sigmab1/K^h3; sigmab2 = sigmab2/K^h1; sigmab3 = sigmab3/K^h2;
para = [h1,h2,h3,d1,d2,d3,rhou1,rhob1,rhou2,rhob2,rhou3,rhob3,...
    sigmau1,sigmab1,sigmau2,sigmab2,sigmau3,sigmab3,M,N,P];

% -------------------------------------------------------------------------
% simulate trajectories
step = 5e-2/K/sigmau1; num = 3e3;
time = 0:step:T; len = length(time);
x1 = zeros(num,len); x2 = zeros(num,len); x3 = zeros(num,len);
s1 = zeros(num,len); s2 = zeros(num,len); s3 = zeros(num,len);
for k = 1:num
    % initial value
    x1(k,1) = 0; x2(k,1) = 0; x3(k,1) = 0;
    s1(k,1) = 1; s2(k,1) = 0; s3(k,1) = 0;
    
    % evolution
    for i = 1:len-1
        rho1 = rhou1*(s1(k,i)==0)+rhob1*(s1(k,i)==1);
        rho2 = rhou2*(s2(k,i)==0)+rhob2*(s2(k,i)==1);
        rho3 = rhou3*(s3(k,i)==0)+rhob3*(s3(k,i)==1);
        rhotot = rho1+rho2+rho3;
        dtot = d1*x1(k,i)+d2*x2(k,i)+d3*x3(k,i);
        prod1 = 1; prod2 = 1; prod3 = 1;
        for ii = 1:h1
            prod1 = prod1*(x1(k,i)-ii+1);
        end
        for ii = 1:h2
            prod2 = prod2*(x2(k,i)-ii+1);
        end
        for ii = 1:h3
            prod3 = prod3*(x3(k,i)-ii+1);
        end
        sigma1 = sigmab1*prod3*(s1(k,i)==0)+sigmau1*(s1(k,i)==1);
        sigma2 = sigmab2*prod1*(s2(k,i)==0)+sigmau2*(s2(k,i)==1);
        sigma3 = sigmab3*prod2*(s3(k,i)==0)+sigmau3*(s3(k,i)==1);
        sigmatot = sigma1+sigma2+sigma3;
        tot = rhotot+dtot+sigmatot;
        if rand < step*tot
            temp = rand;
            if temp < rho1/tot
                x1(k,i+1) = x1(k,i)+1; x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
            elseif temp < (rho1+rho2)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)+1; x3(k,i+1) = x3(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
            elseif temp < rhotot/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i)+1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
            elseif temp < (rhotot+d1*x1(k,i))/tot
                x1(k,i+1) = x1(k,i)-1; x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
            elseif temp < (rhotot+d1*x1(k,i)+d2*x2(k,i))/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-1; x3(k,i+1) = x3(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
            elseif temp < (rhotot+dtot)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i)-1;
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
            elseif temp < (rhotot+dtot+sigma1)/tot
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i)-h3*(s1(k,i)==0)+h3*(s1(k,i)==1);
                s1(k,i+1) = 1-s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
            elseif temp < (rhotot+dtot+sigma1+sigma2)/tot
                x1(k,i+1) = x1(k,i)-h1*(s2(k,i)==0)+h1*(s2(k,i)==1); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = 1-s2(k,i); s3(k,i+1) = s3(k,i);
            else
                x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i)-h2*(s3(k,i)==0)+h2*(s3(k,i)==1); x3(k,i+1) = x3(k,i);
                s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = 1-s3(k,i);
            end
        else
            x1(k,i+1) = x1(k,i); x2(k,i+1) = x2(k,i); x3(k,i+1) = x3(k,i);
            s1(k,i+1) = s1(k,i); s2(k,i+1) = s2(k,i); s3(k,i+1) = s3(k,i);
        end
    end
end

% -------------------------------------------------------------------------
% downsample
numpoint = 200; itv = floor(T/step/numpoint);
time = time(1:itv:end)'; len = length(time);
x1 = x1(:,1:itv:end); x2 = x2(:,1:itv:end); x3 = x3(:,1:itv:end);
s1 = s1(:,1:itv:end); s2 = s2(:,1:itv:end); s3 = s3(:,1:itv:end);

% -------------------------------------------------------------------------
% calculate moments
g0 = sum((s1(:,:)==0))'/num;
g1 = sum((s1(:,:)==1))'/num;
g11 = sum((s1(:,:)==1).*(s2(:,:)==1))'/num;
m0 = sum(x1.*(s1(:,:)==0))'/num;
m1 = sum(x1.*(s1(:,:)==1))'/num;
m10 = sum(x1.*(s1(:,:)==1).*(s2(:,:)==0))'/num;
mm10 = sum(x1.*(x1-1).*(s1(:,:)==1).*(s2(:,:)==0))'/num;
mmm10 = sum(x1.*(x1-1).*(x1-2).*(s1(:,:)==1).*(s2(:,:)==0))'/num;
p0 = sum(x3.*(s1(:,:)==0))'/num;
pp0 = sum(x3.*(x3-1).*(s1(:,:)==0))'/num;
ppp0 = sum(x3.*(x3-1).*(x3-2).*(s1(:,:)==0))'/num;
mp0 = sum(x1.*x3.*(s1(:,:)==0))'/num;
mpp0 = sum(x1.*x3.*(x3-1).*(s1(:,:)==0))'/num;
mppp0 = sum(x1.*x3.*(x3-1).*(x3-2).*(s1(:,:)==0))'/num;
mm0 = sum(x1.*(x1-1).*(s1(:,:)==0))'/num;
mm1 = sum(x1.*(x1-1).*(s1(:,:)==1))'/num;
input = [g0,g1,g11,m0,m1,m10,mm10,mmm10,p0,pp0,ppp0,mp0,mpp0,mppp0,time];

% -------------------------------------------------------------------------
% simulate distributions using SSA
backnum = 0;
temp = abs(time-Tc);
ind = find(temp==min(temp)); ind = min(ind);
distm = zeros(1,M);
for i = 1:M
    distm(i) = sum(sum(x1(:,ind-backnum:ind+backnum)==i-1))/num/(2*backnum+1);
end
hori = 0:maxm-1; vert = distm(1:maxm); dist = distm;
figure; plot(hori,vert,'color',0.6*ones(3,1)); box off; hold on
xmax = 120; xlim([0,xmax]); set(gca,'xtick',[0,xmax/3,xmax*2/3,xmax]);
ymax = 0.033; ylim([0,ymax]); set(gca,'ytick',[0,ymax/3,ymax*2/3,ymax]);
patch([hori fliplr(hori)], [vert 1e-4*ones(size(vert))],[0.9 0.9 0.9],'LineStyle','none'); hold on

% -------------------------------------------------------------------------
% predict distributions using Holimap
birthlma = zeros(1,2*M); birthlma(M+1) = 1;

% distribution predicted by SSA+LMA
sollma = master_LMA(para,input,birthlma,Tc);
distlma = sollma(end,1:M)+sollma(end,M+1:2*M);
distlma = max(distlma,zeros(1,M));
plot(0:maxm-1,distlma,'r','Linewidth',1.5); hold on

% distribution predicted by SSA+2HM
[~,solhm] = ode45(@master_2HM,[0,Tc],birthlma,[],para,input);
disthm = solhm(end,1:M)+solhm(end,M+1:2*M);
disthm = max(disthm,zeros(1,M));
plot(0:maxm-1,disthm(1:maxm),'b','Linewidth',2); hold on

toc