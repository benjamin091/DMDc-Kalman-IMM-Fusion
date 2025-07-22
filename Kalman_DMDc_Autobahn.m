function kal_A = Kalman_DMDc_Autobahn(local_A,n,N_pred,input)
states = 2;
%Systemmatrizen aus DMDc
w = [randn(states*n,length(local_A.Q))]*sqrt(local_A.Q);
v = [randn(states*n,length(local_A.R))]*sqrt(local_A.R);
B = [local_A.res.B,w];
C = [eye(states), zeros(states,(n-1)*states)];
D = 0;

sys = ss(local_A.res.A,B,C,D,local_A.T);


%Bestimmung des Kalman-Modells
[kalmf_Autobahn,L_A,P_A] = kalman(sys,local_A.Q,local_A.R,local_A.N,'current');

sim_A = struct;
t = [1:N_pred];

% sim_A.v = zeros(length(local_A.res.x_pred),N_pred-n+1);
% sim_A.P = zeros(length(local_A.res.x2_pred),N_pred-n+1);
sim_A.v = zeros(length(local_A.res.x_pred),N_pred);
sim_A.P = zeros(length(local_A.res.x2_pred),N_pred);
sim_A.vplot = zeros(length(local_A.res.x_pred),N_pred-n+1);
sim_A.Pplot = zeros(length(local_A.res.x2_pred),N_pred-n+1);

for i = 1:length(local_A.res.x_pred)
    if input == 4
        u = [...
            local_A.res.Y(1,i:N_pred-1+i);...
            local_A.res.Y(2,i:N_pred-1+i);...
            local_A.res.Y(3,i:N_pred-1+i);...
            local_A.res.Y(4,i:N_pred-1+i);...
            local_A.res.x_pred(i,:);...
            local_A.res.x2_pred(i,:)...
            ];
    elseif input == 3
        u = [...
            local_A.res.Y(1,i:N_pred-1+i);...
            local_A.res.Y(2,i:N_pred-1+i);...
            local_A.res.Y(3,i:N_pred-1+i);...
            local_A.res.x_pred(i,:);...
            local_A.res.x2_pred(i,:)...
            ];
    elseif input ==2
        u = [...
            local_A.res.Y(1,i:N_pred-1+i);...
            local_A.res.Y(2,i:N_pred-1+i);...
            local_A.res.x_pred(i,:);...
            local_A.res.x2_pred(i,:)...
            ];
    elseif input == 1
        u = [...
            local_A.res.Y(1,i:N_pred-1+i);...
            local_A.res.x_pred(i,:);...
            local_A.res.x2_pred(i,:)...
            ];
    end
    out = lsim(kalmf_Autobahn,u,t,local_A.res.X(1:states*n,i))';
    
%     sim_A.v(i,:) = out(1,n:end);
%     sim_A.P(i,:) = out(2,n:end);
    sim_A.v(i,1:end) = out(1,:);
    sim_A.P(i,1:end) = out(2,:);
    sim_A.vplot(i,1:end) = out(1,n:end);
    sim_A.Pplot(i,1:end) = out(2,n:end);
%     sim_A.x_A(n*states*(i-1)+1:i*n*states,:) = out(states+1:states*(n+1),n:end);
    sim_A.x_A(n*states*(i-1)+1:i*n*states,:) = out(states+1:states*(n+1),:);
end
    kal_A = struct;
    kal_A.kalmf = kalmf_Autobahn;
    kal_A.P = P_A;
    kal_A.sim = sim_A;
    kal_A.C = C;
end


