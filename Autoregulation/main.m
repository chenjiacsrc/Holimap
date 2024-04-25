warning off all
tic

% -------------------------------------------------------------------------
% model parameters
% non-cooperative binding
% h = 1; d = 1; B = 0.1; rhou = 0.5*d/B; rhob = 20*d/B;
% sigmau = 1; sigmab = 0.5*sigmau;

% cooperative binding
h = 2; d = 1; B = 0.1; rhou = 5*d/B; rhob = 20*d/B;
sigmau = 0.6; sigmab = 0.02*sigmau;

% -------------------------------------------------------------------------
% true distribution computed using FSP
N = 61; maxm = 41; T = 10; Tc = T; len = T*10;
para = [h,d,B,rhou,rhob,sigmau,sigmab,N];
birth = zeros(1,2*N); birth(1) = 1;
time = linspace(0,T,len); itv = T/(len-1);
sol = zeros(len,2*N);
Q = master(para);
R = expm(itv*Q);
sol(1,:) = birth;
for i = 1:len-1
    sol(i+1,:) = sol(i,:)*R;
end
leng = length(para); temp = zeros(len+1,2*N+1);
temp(1:len,1:2*N) = sol; temp(1:len,end) = time';
temp(end,1:leng) = para; para = temp;
temp = abs(time-Tc);
ind = find(temp==min(temp)); ind = min(ind);
distm = sol(ind,1:N)+sol(ind,N+1:2*N);
distm = max(distm,zeros(1,N));
hori = 0:maxm-1; vert = distm(1:maxm);
plot(hori,vert,'color',0.6*ones(3,1)); box off; hold on
xmax = 40; xlim([0,xmax]); set(gca,'xtick',[0,xmax/4,xmax/2,xmax*3/4,xmax]);
patch([hori fliplr(hori)], [vert 1.5e-4*ones(size(vert))],[0.9 0.9 0.9],'LineStyle','none'); hold on

% -------------------------------------------------------------------------
% initial values of moments
momb = zeros(1,7);
momb(1) = sum(birth(1:N));
momb(2) = sum((0:N-1).*birth(1:N));
momb(3) = sum((0:N-1).*birth(N+1:2*N));
momb(4) = sum((0:N-1).*(-1:N-2).*birth(1:N));
momb(5) = sum((0:N-1).*(-1:N-2).*birth(N+1:2*N));
momb(6) = sum((0:N-1).*(-1:N-2).*(-2:N-3).*birth(1:N));
momb(7) = sum((0:N-1).*(-1:N-2).*(-2:N-3).*birth(N+1:2*N));

% -------------------------------------------------------------------------
% distribution predicted by LMA
[time,sol] = ode23(@mom_LMA,[0,Tc],momb(1:5),[],para);
len = length(time);
input = zeros(6,len);
input(1,:) = time; input(2,:) = sol(:,1);
input(3,:) = sol(:,2); input(4,:) = sol(:,3);
input(5,:) = sol(:,4); input(6,:) = sol(:,5);

[~,sol] = ode23(@master_LMA,[0,Tc],birth,[],para,input,Tc);
distlma = sol(end,1:N)+sol(end,N+1:2*N);
distlma = max(distlma,zeros(1,N));
plot(0:maxm-1,distlma(1:maxm),'g','color',[0.47,0.67,0.19]); hold on

% -------------------------------------------------------------------------
% distribution predicted by 2-HM
[time,sol] = ode23(@mom_2HM,[0,Tc],momb,[],para);
len = length(time);
input = zeros(6,len);
input(1,:) = time; input(2,:) = sol(:,1);
input(3,:) = sol(:,2); input(4,:) = sol(:,3);
input(5,:) = sol(:,4); input(6,:) = sol(:,5);
input(7,:) = sol(:,6); input(8,:) = sol(:,7);

[~,sol] = ode23(@master_2HM,[0,Tc],birth,[],para,input);
distlma = sol(end,1:N)+sol(end,N+1:2*N);
distlma = max(distlma,zeros(1,N));
plot(0:maxm-1,distlma(1:maxm),'r'); hold on

% -------------------------------------------------------------------------
% distribution predicted by 2-HM
[time,sol] = ode23(@mom_4HM,[0,Tc],momb,[],para);
len = length(time);
input = zeros(8,len);
input(1,:) = time; input(2,:) = sol(:,1);
input(3,:) = sol(:,2); input(4,:) = sol(:,3);
input(5,:) = sol(:,4); input(6,:) = sol(:,5);
input(7,:) = sol(:,6); input(8,:) = sol(:,7);

[~,sol] = ode23(@master_4HM,[0,Tc],birth,[],para,input);
distlma = sol(end,1:N)+sol(end,N+1:2*N);
distlma = max(distlma,zeros(1,N));
plot(0:maxm-1,distlma(1:maxm),'b'); hold on
legend('FSP','FSP','LMA', '2-HM','4-HM');
legend('boxoff');

toc