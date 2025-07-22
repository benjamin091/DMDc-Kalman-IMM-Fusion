function [IMM,P_IMM] = IMM_DMDc(n,N_pred,markov,mu_hat,local_WC,local_L,local_A,kal_WC,kal_L,kal_A,data)
states = 2;
%% Interaction step
if markov.m == 1
    p = [markov.p11 (1-markov.p11) 0; ...
        ((1-markov.p22)/2) markov.p22 ((1-markov.p22)/2); ...
        0 (1-markov.p33) markov.p33];
elseif markov.m == 2
    p = [markov.p11 (1-markov.p11)/2 (1-markov.p11)/2; ...
        ((1-markov.p22)/2) markov.p22 ((1-markov.p22)/2); ...
        (1-markov.p33)/2 (1-markov.p33)/2 markov.p33];
elseif markov.m == 3
     p = [markov.p11 (1-markov.p11)/2 (1-markov.p11)/2; ...
         0 markov.p22 (1-markov.p22); ...
         0 0 markov.p33];
elseif markov.m == 4
    p = [markov.p11 (1-markov.p11) 0; ...
         0 markov.p22 (1-markov.p22); ...
         0 0 markov.p33];
end
    
    
c_dash = [(mu_hat(1,1)*(p(1,1)+p(2,1)+p(3,1)));...
    (mu_hat(2,1)*(p(1,2)+p(2,2)+p(3,2)));...
    (mu_hat(3,1)*(p(1,3)+p(2,3)+p(3,3)))];

mu_tilde = [p(1,1)*mu_hat(1,1)/c_dash(1,1) p(2,1)*mu_hat(2,1)/c_dash(1,1) p(3,1)*mu_hat(2,1)/c_dash(1,1); ...
    p(1,2)*mu_hat(1,1)/c_dash(2,1) p(2,2)*mu_hat(2,1)/c_dash(2,1) p(3,2)*mu_hat(3,1)/c_dash(2,1); ...
    p(1,3)*mu_hat(1,1)/c_dash(3,1) p(2,3)*mu_hat(3,1)/c_dash(3,1) p(3,3)*mu_hat(3,1)/c_dash(3,1)];

%% blening the state outputs from the last step accourding to
%%the probability of their current contribution to this state
x_tilde_WC = zeros(size(kal_WC.sim.x_WC));
x_tilde_L = zeros(size(kal_L.sim.x_L));
x_tilde_A = zeros(size(kal_A.sim.x_A));

Z_WC = zeros(size(kal_WC.sim.x_WC*states));
Z_L = zeros(size(kal_L.sim.x_L*states));
Z_A = zeros(size(kal_A.sim.x_A*states));

% Lambda_WC = zeros(size(kal_WC.sim.x_WC));
% Lambda_L = zeros(size(kal_L.sim.x_L));
% Lambda_A = zeros(size(kal_A.sim.x_A));

Lambda_WC = zeros(size(kal_WC.sim.x_WC,1)/states/n,size(kal_WC.sim.x_WC,2));
Lambda_L = zeros(size(kal_L.sim.x_L,1)/states/n,size(kal_L.sim.x_L,2));
Lambda_A = zeros(size(kal_A.sim.x_A,1)/states/n,size(kal_A.sim.x_A,2));

mu_IMM_WC = zeros(size(kal_WC.sim.x_WC));
mu_IMM_L = zeros(size(kal_L.sim.x_L));
mu_IMM_A = zeros(size(kal_A.sim.x_A));

x_IMM = zeros(size(kal_WC.sim.x_WC));
x_IMM_Out = zeros(size(kal_WC.sim.x_WC,1)/n,size(kal_WC.sim.x_WC,2));

IMM = struct;
% IMM.v = zeros(size(kal_WC.sim.v));
% IMM.P = zeros(size(kal_WC.sim.v));
IMM.v = zeros(size(kal_WC.sim.vplot));
IMM.P = zeros(size(kal_WC.sim.Pplot));
%% 
for i = 1:length(kal_WC.sim.x_WC)/(states*n)
%     for j = 1:(N_pred-n+1)
    for j = 1:N_pred
    x_tilde_WC(n*states*(i-1)+1:i*n*states,j) = kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)*mu_tilde(1,1)+...
        +kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)*mu_tilde(1,2)+...
        +kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)*mu_tilde(1,3);
    x_tilde_L(n*states*(i-1)+1:i*n*states,j) = kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)*mu_tilde(2,1)+...
        +kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)*mu_tilde(2,2)+...
        +kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)*mu_tilde(2,3);
    x_tilde_A(n*states*(i-1)+1:i*n*states,j) = kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)*mu_tilde(3,1)+...
        +kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)*mu_tilde(3,2)+...
        +kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)*mu_tilde(3,3);

    P_tilde_WC(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states) = mu_tilde(1,1)*(kal_WC.P+(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_tilde_WC(n*states*(i-1)+1:i*n*states,j))*(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_tilde_WC(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_tilde(2,1)*(kal_L.P+(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_tilde_WC(n*states*(i-1)+1:i*n*states,j))*(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_tilde_WC(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_tilde(3,1)*(kal_A.P+(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_tilde_WC(n*states*(i-1)+1:i*n*states,j))*(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_tilde_WC(n*states*(i-1)+1:i*n*states,j))');

    P_tilde_L(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states) = mu_tilde(1,2)*(kal_WC.P+(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_tilde_L(n*states*(i-1)+1:i*n*states,j))*(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_tilde_L(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_tilde(2,2)*(kal_L.P+(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_tilde_L(n*states*(i-1)+1:i*n*states,j))*(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_tilde_L(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_tilde(3,2)*(kal_A.P+(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_tilde_L(n*states*(i-1)+1:i*n*states,j))*(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_tilde_L(n*states*(i-1)+1:i*n*states,j))');

    P_tilde_A(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states) = mu_tilde(1,3)*(kal_WC.P+(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_tilde_A(n*states*(i-1)+1:i*n*states,j))*(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_tilde_A(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_tilde(2,3)*(kal_L.P+(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_tilde_A(n*states*(i-1)+1:i*n*states,j))*(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_tilde_A(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_tilde(3,3)*(kal_A.P+(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_tilde_A(n*states*(i-1)+1:i*n*states,j))*(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_tilde_A(n*states*(i-1)+1:i*n*states,j))');

    Z_WC(states*(i-1)+1,j) = data.v(i,j) - kal_WC.sim.v(i,j);
    Z_WC(states*(i-1)+2,j) = data.P(i,j) - kal_WC.sim.P(i,j);
%     Z_WC(states*i-1,j) = data.adot(i,j) - kal_WC.sim.adot(i,j);
%     Z_WC(states*i,j) = data.a(i,j) - kal_WC.sim.a(i,j);
    
    Z_L(states*(i-1)+1,j) = data.v(i,j) - kal_L.sim.v(i,j);
    Z_L(states*(i-1)+2,j) = data.P(i,j) - kal_L.sim.P(i,j);
%     Z_L(states*i-1,j) = data.adot(i,j) - kal_L.sim.adot(i,j);
%     Z_L(states*i,j) = data.a(i,j) - kal_L.sim.a(i,j);
    
    Z_A(states*(i-1)+1,j) = data.v(i,j) - kal_A.sim.v(i,j);
    Z_A(states*(i-1)+2,j) = data.P(i,j) - kal_A.sim.P(i,j);
%     Z_A(states*i-1,j) = data.adot(i,j) - kal_A.sim.adot(i,j);
%     Z_A(states*i,j) = data.a(i,j) - kal_A.sim.a(i,j);

    S_WC(states*(i-1)+1:i*states,states*(j-1)+1:j*states) = kal_WC.C*P_tilde_WC(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states)*kal_WC.C'+local_WC.R;
    S_L(states*(i-1)+1:i*states,states*(j-1)+1:j*states) = kal_L.C*P_tilde_L(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states)*kal_L.C'+local_L.R;
    S_A(states*(i-1)+1:i*states,states*(j-1)+1:j*states) = kal_A.C*P_tilde_A(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states)*kal_A.C'+local_A.R;

    Lambda_WC(i,j) = sqrt(det(2*pi*S_WC(states*(i-1)+1:i*states,states*(j-1)+1:j*states)))^-1*exp(-.5*Z_WC(states*(i-1)+1:i*states,j)'/S_WC(states*(i-1)+1:i*states,states*(j-1)+1:j*states)*Z_WC(states*(i-1)+1:i*states,j));
    Lambda_L(i,j) = sqrt(det(2*pi*S_L(states*(i-1)+1:i*states,states*(j-1)+1:j*states)))^-1*exp(-.5*Z_L(states*(i-1)+1:i*states,j)'/S_L(states*(i-1)+1:i*states,states*(j-1)+1:j*states)*Z_L(states*(i-1)+1:i*states,j));
    Lambda_A(i,j) = sqrt(det(2*pi*S_A(states*(i-1)+1:i*states,states*(j-1)+1:j*states)))^-1*exp(-.5*Z_A(states*(i-1)+1:i*states,j)'/S_A(states*(i-1)+1:i*states,states*(j-1)+1:j*states)*Z_A(states*(i-1)+1:i*states,j));

    mu_IMM_WC(i,j) = Lambda_WC(i,j)*c_dash(1,1)/(Lambda_WC(i,j)*c_dash(1,1)+Lambda_L(i,j)*c_dash(2,1)+Lambda_A(i,j)*c_dash(3,1));
    mu_IMM_L(i,j) = Lambda_L(i,j)*c_dash(2,1)/(Lambda_WC(i,j)*c_dash(1,1)+Lambda_L(i,j)*c_dash(2,1)+Lambda_A(i,j)*c_dash(3,1));
    mu_IMM_A(i,j) = Lambda_A(i,j)*c_dash(3,1)/(Lambda_WC(i,j)*c_dash(1,1)+Lambda_L(i,j)*c_dash(2,1)+Lambda_A(i,j)*c_dash(3,1));

    x_IMM(n*states*(i-1)+1:i*n*states,j) = kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)*mu_IMM_WC(i,j)+kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)*mu_IMM_L(i,j)+kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)*mu_IMM_A(i,j);
    
    P_IMM(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states) = mu_IMM_WC(i,j)*(P_tilde_WC(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states)+(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_IMM(n*states*(i-1)+1:i*n*states,j))*(kal_WC.sim.x_WC(n*states*(i-1)+1:i*n*states,j)-x_IMM(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_IMM_L(i,j)*(P_tilde_L(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states)+(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_IMM(n*states*(i-1)+1:i*n*states,j))*(kal_L.sim.x_L(n*states*(i-1)+1:i*n*states,j)-x_IMM(n*states*(i-1)+1:i*n*states,j))') + ...
        + mu_IMM_A(i,j)*(P_tilde_A(n*states*(i-1)+1:i*n*states,n*states*(j-1)+1:j*n*states)+(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_IMM(n*states*(i-1)+1:i*n*states,j))*(kal_A.sim.x_A(n*states*(i-1)+1:i*n*states,j)-x_IMM(n*states*(i-1)+1:i*n*states,j))');
    end
end

for i = 1:(length(kal_WC.sim.v))
%     for j = 1:(N_pred-n+1)
    for j = 1:(N_pred)
        x_IMM_Out(states*(i-1)+1:states*i,j) = x_IMM(states*n*(i-1)+1:states*n*(i-1)+states,j);
    end
%     IMM.v(i,:) = x_IMM_Out(states*(i-1)+1,:);
%     IMM.P(i,:) = x_IMM_Out(states*(i-1)+2,:);
    IMM.v(i,:) = x_IMM_Out(states*(i-1)+1,n:N_pred);
    IMM.P(i,:) = x_IMM_Out(states*(i-1)+2,n:N_pred);
end
end