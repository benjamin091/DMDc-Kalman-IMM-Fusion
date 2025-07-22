function kal_L = Kalman_DMDc_Landstrasse(local_L,n,N_pred,input)
states = 2;
%Systemmatrizen aus DMDc
w = [randn(states*n,length(local_L.Q))]*sqrt(local_L.Q);
v = [randn(states*n,length(local_L.R))]*sqrt(local_L.R);
B = [local_L.res.B,w];
C = [eye(states), zeros(states,(n-1)*states)];
D = 0;

sys = ss(local_L.res.A,B,C,D,local_L.T);


%Bestimmung des Kalman-Modells
[kalmf_Landstrasse,L_L,P_L] = kalman(sys,local_L.Q,local_L.R,local_L.N,'current');

sim_L = struct;
t = [1:N_pred];

% sim_L.v = zeros(length(local_L.res.x_pred),N_pred-n+1);
% sim_L.P = zeros(length(local_L.res.x2_pred),N_pred-n+1);
sim_L.v = zeros(length(local_L.res.x_pred),N_pred);
    sim_L.P = zeros(length(local_L.res.x2_pred),N_pred);
    sim_L.vplot = zeros(length(local_L.res.x_pred),N_pred-n+1);
    sim_L.Pplot = zeros(length(local_L.res.x2_pred),N_pred-n+1);
for i = 1:length(local_L.res.x_pred)
    if input ==4
        u = [...
            local_L.res.Y(1,i:N_pred-1+i);...
            local_L.res.Y(2,i:N_pred-1+i);...
            local_L.res.Y(3,i:N_pred-1+i);...
            local_L.res.Y(4,i:N_pred-1+i);...
            local_L.res.x_pred(i,:);...
            local_L.res.x2_pred(i,:)...
            ];
    elseif input == 3
        u = [...
            local_L.res.Y(1,i:N_pred-1+i);...
            local_L.res.Y(2,i:N_pred-1+i);...
            local_L.res.Y(3,i:N_pred-1+i);...
            local_L.res.x_pred(i,:);...
            local_L.res.x2_pred(i,:)...
            ];
    elseif input ==2
        u = [...
            local_L.res.Y(1,i:N_pred-1+i);...
            local_L.res.Y(2,i:N_pred-1+i);...
            local_L.res.x_pred(i,:);...
            local_L.res.x2_pred(i,:)...
            ];
    elseif input == 1
        u = [...
            local_L.res.Y(1,i:N_pred-1+i);...
            local_L.res.x_pred(i,:);...
            local_L.res.x2_pred(i,:)...
            ];
    end
    out = lsim(kalmf_Landstrasse,u,t,local_L.res.X(1:states*n,i))';
    
%     sim_L.v(i,:) = out(1,n:end);
%     sim_L.P(i,:) = out(2,n:end);
    sim_L.v(i,1:end) = out(1,:);
    sim_L.P(i,1:end) = out(2,:);
    sim_L.vplot(i,1:end) = out(1,n:end);
    sim_L.Pplot(i,1:end) = out(2,n:end);
%     sim_L.x_L(n*states*(i-1)+1:i*n*states,:) = out(states+1:states*(n+1),n:end);
    sim_L.x_L(n*states*(i-1)+1:i*n*states,:) = out(states+1:states*(n+1),:);
end
    kal_L = struct;
    kal_L.kalmf = kalmf_Landstrasse;
    kal_L.P = P_L;
    kal_L.sim = sim_L;
    kal_L.C = C;
end


