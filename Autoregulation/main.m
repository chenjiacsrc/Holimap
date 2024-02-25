warning off all
tic

% -------------------------------------------------------------------------
% model parameters
h = 2; d = 1; B = 0.1; rhou = 0.5*d/B; rhob = 20*d/B;
sigmau = 1; sigmab = 1;

% -------------------------------------------------------------------------
% true distribution computed using FSP
N = 61; maxm = 41; T = 10; Tc = 10; len = T*10;
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

inputg0 = zeros(1,len);
inputp0 = zeros(1,len); inputp1 = zeros(1,len);
inputpp0 = zeros(1,len); inputpp1 = zeros(1,len);
inputppp0 = zeros(1,len); inputppp1 = zeros(1,len);
for k = 1:len
    inputg0(k) = sum(sol(k,1:N));
    for i = 0:N-1
        inputp0(k) = inputp0(k)+i*sol(k,i+1);
        inputp1(k) = inputp1(k)+i*sol(k,i+1+N);
        inputpp0(k) = inputpp0(k)+i*(i-1)*sol(k,i+1);
        inputpp1(k) = inputpp1(k)+i*(i-1)*sol(k,i+1+N);
        inputppp0(k) = inputppp0(k)+i*(i-1)*(i-2)*sol(k,i+1);
        inputppp1(k) = inputppp1(k)+i*(i-1)*(i-2)*sol(k,i+1+N);
    end
end
input = zeros(6,len);
input(1,:) = time; input(2,:) = inputg0;
input(3,:) = inputp0; input(4,:) = inputp1;
input(5,:) = inputpp0; input(6,:) = inputpp1;
input(7,:) = inputppp0; input(8,:) = inputppp1;

temp = abs(time-Tc);
ind = find(temp==min(temp)); ind = min(ind);
distm = sol(ind,1:N)+sol(ind,N+1:2*N);
distm = max(distm,zeros(1,N));
hori = 0:maxm-1; vert = distm(1:maxm);
figure; plot(hori,vert,'color',0.6*ones(3,1));
box on; hold on
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
inputt = zeros(6,len);
inputt(1,:) = time; inputt(2,:) = sol(:,1);
inputt(3,:) = sol(:,2); inputt(4,:) = sol(:,3);
inputt(5,:) = sol(:,4); inputt(6,:) = sol(:,5);

[~,sol] = ode23(@master_LMA,[0,Tc],birth,[],para,inputt,Tc);
distlma = sol(end,1:N)+sol(end,N+1:2*N);
distlma = max(distlma,zeros(1,N));
plot(0:maxm-1,distlma(1:maxm),'r'); hold on

% -------------------------------------------------------------------------
% distribution predicted by 2-HM
[time,sol] = ode23(@mom_2HM,[0,Tc],momb,[],para);
len = length(time);
inputt = zeros(6,len);
inputt(1,:) = time; inputt(2,:) = sol(:,1);
inputt(3,:) = sol(:,2); inputt(4,:) = sol(:,3);
inputt(5,:) = sol(:,4); inputt(6,:) = sol(:,5);
inputt(7,:) = sol(:,6); inputt(8,:) = sol(:,7);

[~,sol] = ode23(@master_2HM,[0,Tc],birth,[],para,input);
distlma = sol(end,1:N)+sol(end,N+1:2*N);
distlma = max(distlma,zeros(1,N));
plot(0:2:maxm-1,distlma(1:2:maxm),'.','MarkerSize',18,'color',[0.47,0.67,0.19]); hold on

% -------------------------------------------------------------------------
% distribution predicted by 4-HM
[time,sol] = ode23(@mom_4HM,[0,Tc],momb,[],para);
len = length(time);
inputt = zeros(6,len);
inputt(1,:) = time; inputt(2,:) = sol(:,1);
inputt(3,:) = sol(:,2); inputt(4,:) = sol(:,3);
inputt(5,:) = sol(:,4); inputt(6,:) = sol(:,5);
inputt(7,:) = sol(:,6); inputt(8,:) = sol(:,7);

[~,sol] = ode23(@master_4HM,[0,Tc],birth,[],para,input);
distlma = sol(end,1:N)+sol(end,N+1:2*N);
distlma = max(distlma,zeros(1,N));
plot(0:maxm-1,distlma(1:maxm),'b'); hold on

toc