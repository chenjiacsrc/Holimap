function f = mom_LMA(t,x,para)

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

% -------------------------------------------------------------------------
% moment dynamics
if h == 1
    sigmabbar = sigmab*p0/g0;
else
    sigmabbar = sigmab*pp0/g0;
end
% dynamics of gene state
f01 = sigmau*g1-sigmabbar*g0;
% dynmaics of first moments
f02 = rhou*B*g0-d*p0+sigmau*p1-sigmabbar*p0;
f03 = rhob*B*g1-d*p1+sigmabbar*p0-sigmau*p1;
% dynmaics of second moments
f04 = 2*rhou*B*(p0+B*g0)-2*d*pp0+sigmau*pp1-sigmabbar*pp0;
f05 = 2*rhob*B*(p1+B*g1)-2*d*pp1+sigmabbar*pp0-sigmau*pp1;
f = [f01,f02,f03,f04,f05]';