function kal_WC = Kalman_DMDc_WienCity(local_WC,n,N_pred,input)
states = 2;
%Systemmatrizen aus DMDc
w = [randn(states*n,length(local_WC.Q))]*sqrt(local_WC.Q);
v = [randn(states*n,length(local_WC.R))]*sqrt(local_WC.R);
B = [local_WC.res.B,w];
C = [eye(states), zeros(states,(n-1)*states)];
D = 0;

sys = ss(local_WC.res.A,B,C,D,local_WC.T);


%Bestimmung des Kalman-Modells
[kalmf_WienCity,L_WC,P_WC] = kalman(sys,local_WC.Q,local_WC.R,local_WC.N,'current');

sim_WC = struct;
t = [1:N_pred];
% sim_WC.v = zeros(length(local_WC.res.x_pred),N_pred-n+1);
% sim_WC.P = zeros(length(local_WC.res.x2_pred),N_pred-n+1);
    sim_WC.v = zeros(length(local_WC.res.x_pred),N_pred);
    sim_WC.P = zeros(length(local_WC.res.x2_pred),N_pred);
    sim_WC.vplot = zeros(length(local_WC.res.x_pred),N_pred-n+1);
    sim_WC.Pplot = zeros(length(local_WC.res.x2_pred),N_pred-n+1);
for i = 1:length(local_WC.res.x_pred)
    if input == 4
        u = [...
            local_WC.res.Y(1,i:N_pred-1+i);...
            local_WC.res.Y(2,i:N_pred-1+i);...
            local_WC.res.Y(3,i:N_pred-1+i);...
            local_WC.res.Y(4,i:N_pred-1+i);...
            local_WC.res.x_pred(i,:); ...
            local_WC.res.x2_pred(i,:)...
            ];
    elseif input == 3
        u = [...
            local_WC.res.Y(1,i:N_pred-1+i);...
            local_WC.res.Y(2,i:N_pred-1+i);...
            local_WC.res.Y(3,i:N_pred-1+i);...
            local_WC.res.x_pred(i,:);...
            local_WC.res.x2_pred(i,:)...
            ];
    elseif input ==2
        u = [...
            local_WC.res.Y(1,i:N_pred-1+i);...
            local_WC.res.Y(2,i:N_pred-1+i);...
            local_WC.res.x_pred(i,:);...
            local_WC.res.x2_pred(i,:)...
            ];
    elseif input == 1
        u = [...
            local_WC.res.Y(1,i:N_pred-1+i);...
            local_WC.res.x_pred(i,:);...
            local_WC.res.x2_pred(i,:)...
            ];
    end
    out = lsim(kalmf_WienCity,u,t,local_WC.res.X(1:states*n,i))';
    
%     sim_WC.v(i,:) = out(1,n:end);
%     sim_WC.P(i,:) = out(2,n:end);
    sim_WC.v(i,1:end) = out(1,:);
    sim_WC.P(i,1:end) = out(2,:);
    sim_WC.vplot(i,1:end) = out(1,n:end);
    sim_WC.Pplot(i,1:end) = out(2,n:end);
%     sim_WC.x_WC(n*states*(i-1)+1:i*n*states,:) = out(states+1:states*(n+1),n:end);
    sim_WC.x_WC(n*states*(i-1)+1:i*n*states,:) = out(states+1:states*(n+1),:);
end
    kal_WC = struct;
    kal_WC.kalmf = kalmf_WienCity;
    kal_WC.P = P_WC;
    kal_WC.sim = sim_WC;
    kal_WC.C = C;
end


