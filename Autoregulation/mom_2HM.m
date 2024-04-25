function f = mom_2HM(t,x,para)

% -------------------------------------------------------------------------
% model parameters
h = para(end,1); d = para(end,2); B = para(end,3);
rhou = para(end,4); rhob = para(end,5);
sigmau = para(end,6); sigmab = para(end,7);                                                                                                             

% -------------------------------------------------------------------------
% variables
eps = 1e-5;
g0 = x(1)+eps; g1 = 1-g0+eps;
p0 = x(2)+eps; p1 = x(3)+eps;
pp0 = x(4); pp1 = x(5)+eps;
ppp0 = x(6)+eps; ppp1 = x(7)+eps;

% ---------------------------------------------------------------------
% calculate the effective binding and unbinding rates
mat = zeros(2,2);
mat(1,:) = [g1,-g0];
mat(2,:) = [p1,-p0];
aim = zeros(2,1);
if h == 1
    aim(1) = sigmau*g1-sigmab*p0;
    aim(2) = sigmau*p1-sigmab*pp0;
else
    aim(1) = sigmau*g1-sigmab*pp0;
    aim(2) = sigmau*p1-sigmab*ppp0;
end
temp = inv(mat)*aim;
sigmaubar = temp(1);
sigmabbar = temp(2);

% -------------------------------------------------------------------------
% moment dynamics
f01 = sigmaubar*g1-sigmabbar*g0;
f02 = rhou*B*g0-d*p0+sigmaubar*p1-sigmabbar*p0;
f03 = rhob*B*g1-d*p1+sigmabbar*p0-sigmaubar*p1;
f04 = 2*rhou*B*(p0+B*g0)-2*d*pp0+sigmaubar*pp1-sigmabbar*pp0;
f05 = 2*rhob*B*(p1+B*g1)-2*d*pp1+sigmabbar*pp0-sigmaubar*pp1;
f06 = 3*rhou*B*(pp0+2*B*p0+2*B^2*g0)-3*d*ppp0+sigmaubar*ppp1-sigmabbar*ppp0;
f07 = 3*rhob*B*(pp1+2*B*p1+2*B^2*g1)-3*d*ppp1+sigmabbar*ppp0-sigmaubar*ppp1;
f = [f01,f02,f03,f04,f05,f06,f07]';
