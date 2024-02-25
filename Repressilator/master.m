function f = master(para)

% -------------------------------------------------------------------------
% model parameters
h1 = para(1); h2 = para(2); h3 = para(3);
d1 = para(4); d2 = para(5); d3 = para(6);
rhou1 = para(7); rhob1 = para(8);
rhou2 = para(9); rhob2 = para(10);
rhou3 = para(11); rhob3 = para(12);
sigmau1 = para(13); sigmab1 = para(14);
sigmau2 = para(15); sigmab2 = para(16); 
sigmau3 = para(17); sigmab3 = para(18);
M = para(19); N = para(20); P = para(21);

% -------------------------------------------------------------------------
% generator matrix
num = 8*M*N*P;
Q = zeros(num);
one = ones(num,1);
% synthesis and degradation the first protein
for i = 0:M-1
    for j = 0:N-1
        for k = 0:P-1
            if i <= M-2
                Q(i*N*P+j*P+k+1+0*M*N*P,(i+1)*N*P+j*P+k+1+0*M*N*P) = rhou1;
                Q(i*N*P+j*P+k+1+1*M*N*P,(i+1)*N*P+j*P+k+1+1*M*N*P) = rhou1;
                Q(i*N*P+j*P+k+1+2*M*N*P,(i+1)*N*P+j*P+k+1+2*M*N*P) = rhou1;
                Q(i*N*P+j*P+k+1+3*M*N*P,(i+1)*N*P+j*P+k+1+3*M*N*P) = rhou1;
                Q(i*N*P+j*P+k+1+4*M*N*P,(i+1)*N*P+j*P+k+1+4*M*N*P) = rhob1;
                Q(i*N*P+j*P+k+1+5*M*N*P,(i+1)*N*P+j*P+k+1+5*M*N*P) = rhob1;
                Q(i*N*P+j*P+k+1+6*M*N*P,(i+1)*N*P+j*P+k+1+6*M*N*P) = rhob1;
                Q(i*N*P+j*P+k+1+7*M*N*P,(i+1)*N*P+j*P+k+1+7*M*N*P) = rhob1;
                Q((i+1)*N*P+j*P+k+1+0*M*N*P,i*N*P+j*P+k+1+0*M*N*P) = (i+1)*d1;
                Q((i+1)*N*P+j*P+k+1+1*M*N*P,i*N*P+j*P+k+1+1*M*N*P) = (i+1)*d1;
                Q((i+1)*N*P+j*P+k+1+2*M*N*P,i*N*P+j*P+k+1+2*M*N*P) = (i+1)*d1;
                Q((i+1)*N*P+j*P+k+1+3*M*N*P,i*N*P+j*P+k+1+3*M*N*P) = (i+1)*d1;
                Q((i+1)*N*P+j*P+k+1+4*M*N*P,i*N*P+j*P+k+1+4*M*N*P) = (i+1)*d1;
                Q((i+1)*N*P+j*P+k+1+5*M*N*P,i*N*P+j*P+k+1+5*M*N*P) = (i+1)*d1;
                Q((i+1)*N*P+j*P+k+1+6*M*N*P,i*N*P+j*P+k+1+6*M*N*P) = (i+1)*d1;
                Q((i+1)*N*P+j*P+k+1+7*M*N*P,i*N*P+j*P+k+1+7*M*N*P) = (i+1)*d1;
            end
            
            % synthesis and degradation the second protein
            if j <= N-2
                Q(i*N*P+j*P+k+1+0*M*N*P,i*N*P+(j+1)*P+k+1+0*M*N*P) = rhou2;
                Q(i*N*P+j*P+k+1+1*M*N*P,i*N*P+(j+1)*P+k+1+1*M*N*P) = rhou2;
                Q(i*N*P+j*P+k+1+2*M*N*P,i*N*P+(j+1)*P+k+1+2*M*N*P) = rhob2;
                Q(i*N*P+j*P+k+1+3*M*N*P,i*N*P+(j+1)*P+k+1+3*M*N*P) = rhob2;
                Q(i*N*P+j*P+k+1+4*M*N*P,i*N*P+(j+1)*P+k+1+4*M*N*P) = rhou2;
                Q(i*N*P+j*P+k+1+5*M*N*P,i*N*P+(j+1)*P+k+1+5*M*N*P) = rhou2;
                Q(i*N*P+j*P+k+1+6*M*N*P,i*N*P+(j+1)*P+k+1+6*M*N*P) = rhob2;
                Q(i*N*P+j*P+k+1+7*M*N*P,i*N*P+(j+1)*P+k+1+7*M*N*P) = rhob2;
                Q(i*N*P+(j+1)*P+k+1+0*M*N*P,i*N*P+j*P+k+1+0*M*N*P) = (j+1)*d2;
                Q(i*N*P+(j+1)*P+k+1+1*M*N*P,i*N*P+j*P+k+1+1*M*N*P) = (j+1)*d2;
                Q(i*N*P+(j+1)*P+k+1+2*M*N*P,i*N*P+j*P+k+1+2*M*N*P) = (j+1)*d2;
                Q(i*N*P+(j+1)*P+k+1+3*M*N*P,i*N*P+j*P+k+1+3*M*N*P) = (j+1)*d2;
                Q(i*N*P+(j+1)*P+k+1+4*M*N*P,i*N*P+j*P+k+1+4*M*N*P) = (j+1)*d2;
                Q(i*N*P+(j+1)*P+k+1+5*M*N*P,i*N*P+j*P+k+1+5*M*N*P) = (j+1)*d2;
                Q(i*N*P+(j+1)*P+k+1+6*M*N*P,i*N*P+j*P+k+1+6*M*N*P) = (j+1)*d2;
                Q(i*N*P+(j+1)*P+k+1+7*M*N*P,i*N*P+j*P+k+1+7*M*N*P) = (j+1)*d2;
            end
            
            % synthesis and degradation the third protein
            if k <= P-2
                Q(i*N*P+j*P+k+1+0*M*N*P,i*N*P+j*P+k+2+0*M*N*P) = rhou3;
                Q(i*N*P+j*P+k+1+1*M*N*P,i*N*P+j*P+k+2+1*M*N*P) = rhob3;
                Q(i*N*P+j*P+k+1+2*M*N*P,i*N*P+j*P+k+2+2*M*N*P) = rhou3;
                Q(i*N*P+j*P+k+1+3*M*N*P,i*N*P+j*P+k+2+3*M*N*P) = rhob3;
                Q(i*N*P+j*P+k+1+4*M*N*P,i*N*P+j*P+k+2+4*M*N*P) = rhou3;
                Q(i*N*P+j*P+k+1+5*M*N*P,i*N*P+j*P+k+2+5*M*N*P) = rhob3;
                Q(i*N*P+j*P+k+1+6*M*N*P,i*N*P+j*P+k+2+6*M*N*P) = rhou3;
                Q(i*N*P+j*P+k+1+7*M*N*P,i*N*P+j*P+k+2+7*M*N*P) = rhob3;
                Q(i*N*P+j*P+k+2+0*M*N*P,i*N*P+j*P+k+1+0*M*N*P) = (k+1)*d3;
                Q(i*N*P+j*P+k+2+1*M*N*P,i*N*P+j*P+k+1+1*M*N*P) = (k+1)*d3;
                Q(i*N*P+j*P+k+2+2*M*N*P,i*N*P+j*P+k+1+2*M*N*P) = (k+1)*d3;
                Q(i*N*P+j*P+k+2+3*M*N*P,i*N*P+j*P+k+1+3*M*N*P) = (k+1)*d3;
                Q(i*N*P+j*P+k+2+4*M*N*P,i*N*P+j*P+k+1+4*M*N*P) = (k+1)*d3;
                Q(i*N*P+j*P+k+2+5*M*N*P,i*N*P+j*P+k+1+5*M*N*P) = (k+1)*d3;
                Q(i*N*P+j*P+k+2+6*M*N*P,i*N*P+j*P+k+1+6*M*N*P) = (k+1)*d3;
                Q(i*N*P+j*P+k+2+7*M*N*P,i*N*P+j*P+k+1+7*M*N*P) = (k+1)*d3;
            end
        end
    end
end
% gene state switching
for i = 0:M-1-h1
    for j = 0:N-1-h2
        for k = 0:P-1-h3
            prop1 = prod(i+1:i+h1);
            prop2 = prod(j+1:j+h2);
            prop3 = prod(k+1:k+h3);
            Q(i*N*P+j*P+k+h3+1+0*M*N*P,i*N*P+j*P+k+1+4*M*N*P) = sigmab1*prop3;
            Q(i*N*P+j*P+k+1+4*M*N*P,i*N*P+j*P+k+h3+1+0*M*N*P) = sigmau1;
            Q(i*N*P+j*P+k+h3+1+1*M*N*P,i*N*P+j*P+k+1+5*M*N*P) = sigmab1*prop3;
            Q(i*N*P+j*P+k+1+5*M*N*P,i*N*P+j*P+k+h3+1+1*M*N*P) = sigmau1;
            Q(i*N*P+j*P+k+h3+1+2*M*N*P,i*N*P+j*P+k+1+6*M*N*P) = sigmab1*prop3;
            Q(i*N*P+j*P+k+1+6*M*N*P,i*N*P+j*P+k+h3+1+2*M*N*P) = sigmau1;
            Q(i*N*P+j*P+k+h3+1+3*M*N*P,i*N*P+j*P+k+1+7*M*N*P) = sigmab1*prop3;
            Q(i*N*P+j*P+k+1+7*M*N*P,i*N*P+j*P+k+h3+1+3*M*N*P) = sigmau1;
            Q((i+h1)*N*P+j*P+k+1+0*M*N*P,i*N*P+j*P+k+1+2*M*N*P) = sigmab2*prop1;
            Q(i*N*P+j*P+k+1+2*M*N*P,(i+h1)*N*P+j*P+k+1+0*M*N*P) = sigmau2;
            Q((i+h1)*N*P+j*P+k+1+1*M*N*P,i*N*P+j*P+k+1+3*M*N*P) = sigmab2*prop1;
            Q(i*N*P+j*P+k+1+3*M*N*P,(i+h1)*N*P+j*P+k+1+1*M*N*P) = sigmau2;
            Q((i+h1)*N*P+j*P+k+1+4*M*N*P,i*N*P+j*P+k+1+6*M*N*P) = sigmab2*prop1;
            Q(i*N*P+j*P+k+1+6*M*N*P,(i+h1)*N*P+j*P+k+1+4*M*N*P) = sigmau2;
            Q((i+h1)*N*P+j*P+k+1+5*M*N*P,i*N*P+j*P+k+1+7*M*N*P) = sigmab2*prop1;
            Q(i*N*P+j*P+k+1+7*M*N*P,(i+h1)*N*P+j*P+k+1+5*M*N*P) = sigmau2;
            Q(i*N*P+(j+h2)*P+k+1+0*M*N*P,i*N*P+j*P+k+1+1*M*N*P) = sigmab3*prop2;
            Q(i*N*P+j*P+k+1+1*M*N*P,i*N*P+(j+h2)*P+k+1+0*M*N*P) = sigmau3;
            Q(i*N*P+(j+h2)*P+k+1+2*M*N*P,i*N*P+j*P+k+1+3*M*N*P) = sigmab3*prop2;
            Q(i*N*P+j*P+k+1+3*M*N*P,i*N*P+(j+h2)*P+k+1+2*M*N*P) = sigmau3;
            Q(i*N*P+(j+h2)*P+k+1+4*M*N*P,i*N*P+j*P+k+1+5*M*N*P) = sigmab3*prop2;
            Q(i*N*P+j*P+k+1+5*M*N*P,i*N*P+(j+h2)*P+k+1+4*M*N*P) = sigmau3;
            Q(i*N*P+(j+h2)*P+k+1+6*M*N*P,i*N*P+j*P+k+1+7*M*N*P) = sigmab3*prop2;
            Q(i*N*P+j*P+k+1+7*M*N*P,i*N*P+(j+h2)*P+k+1+6*M*N*P) = sigmau3;
        end
    end
end
temp = Q*one;
for i = 1:num
    Q(i,i) = -temp(i);
end
f = Q;