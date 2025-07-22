function data = data_fct(input_WC,input_L,input_A,local_WC,local_L,local_A,N_pred,n,states,input,data_WC,data_L,data_A)
data = struct;
%%%Anteile an Strassentypen
data.WC = data_WC/3;
data.L = data_L/3;
data.A = data_A/3;
%%%Vorbereitung der Innovationen für das IMM
% data.v = zeros(length(local_WC.res.x_pred),N_pred-n+1);
% data.P = zeros(length(local_WC.res.x2_pred),N_pred-n+1);
% data.v1 = zeros(length(local_WC.res.x_pred),N_pred-n+1);
% data.P1 = zeros(length(local_WC.res.x2_pred),N_pred-n+1);
% data.v2 = zeros(length(local_WC.res.x_pred),N_pred-n+1);
% data.P2 = zeros(length(local_WC.res.x2_pred),N_pred-n+1);
% data.v3 = zeros(length(local_WC.res.x_pred),N_pred-n+1);
% data.P3 = zeros(length(local_WC.res.x2_pred),N_pred-n+1);
data.v = zeros(size(local_WC.res.x_pred));
data.P = zeros(size(local_WC.res.x_pred));
data.v1 = zeros(size(local_WC.res.x_pred));
data.P1 = zeros(size(local_WC.res.x_pred));
data.v2 = zeros(size(local_WC.res.x_pred));
data.P2 = zeros(size(local_WC.res.x_pred));
data.v3 = zeros(size(local_WC.res.x_pred));
data.P3 = zeros(size(local_WC.res.x_pred));


for i = 1:(length(input_WC.v))-n*states-n-1-N_pred
%         data.v1(i,1:N_pred-n+1) = input_WC.v(n-1+i:i+N_pred-1,:);
%         data.P1(i,1:N_pred-n+1) = input_WC.P(n-1+i:i+N_pred-1,:);
%         
%         data.v2(i,1:N_pred-n+1) = input_L.v(n-1+i:i+N_pred-1,:);
%         data.P2(i,1:N_pred-n+1) = input_L.P(n-1+i:i+N_pred-1,:);
%         
%         data.v3(i,1:N_pred-n+1) = input_A.v(n-1+i:i+N_pred-1,:);
%         data.P3(i,1:N_pred-n+1) = input_A.P(n-1+i:i+N_pred-1,:);
%         

        data.v1(i,1:N_pred) = input_WC.v(i:i+N_pred-1,:);
        data.P1(i,1:N_pred) = input_WC.P(i:i+N_pred-1,:);
        
        data.v2(i,1:N_pred) = input_L.v(i:i+N_pred-1,:);
        data.P2(i,1:N_pred) = input_L.P(i:i+N_pred-1,:);
        
        data.v3(i,1:N_pred) = input_A.v(i:i+N_pred-1,:);
        data.P3(i,1:N_pred) = input_A.P(i:i+N_pred-1,:);
        
        if states == 3
            data.a1(i,1:N_pred) = input_WC.a(i:i+N_pred-1,:);
            data.a2(i,1:N_pred) = input_L.a(i:i+N_pred-1,:);
            data.a3(i,1:N_pred) = input_A.a(i:i+N_pred-1,:);
        elseif states == 4
            data.a1(i,1:N_pred) = input_WC.a(i:i+N_pred-1,:);
            data.a2(i,1:N_pred) = input_L.a(i:i+N_pred-1,:);
            data.a3(i,1:N_pred) = input_A.a(i:i+N_pred-1,:);
            
            data.adot1(i,1:N_pred) = input_WC.adot(i:i+N_pred-1,:);
            data.adot2(i,1:N_pred) = input_L.adot(i:i+N_pred-1,:);
            data.adot3(i,1:N_pred) = input_A.adot(i:i+N_pred-1,:);
        end
end

    data.v = cat(1,data.v1(1:length(data.v1)*data.WC,:),data.v2(length(data.v2)*data.WC+1:length(data.v2)*(data.L+data.WC),:),data.v3(length(data.v3)*(data.L+data.WC)+1:length(data.v3),:));
    data.P = cat(1,data.P1(1:length(data.P1)*data.WC,:),data.P2(length(data.P2)*data.WC+1:length(data.P2)*(data.L+data.WC),:),data.P3(length(data.P3)*(data.L+data.WC)+1:length(data.P3),:));
    
%%%Vorbereitung für die Referenzdaten für den Plot
%%Daten
    data.vref = cat(1,input_WC.v(1:length(local_WC.res.x_pred)*data.WC),input_L.v(length(local_WC.res.x_pred)*data.WC+1:length(local_WC.res.x_pred)*(data.WC+data.L)),input_A.v(length(local_WC.res.x_pred)*(data.WC+data.L)+1:length(local_WC.res.x_pred)));
    data.Pref = cat(1,input_WC.P(1:length(local_WC.res.x_pred)*data.WC),input_L.P(length(local_WC.res.x_pred)*data.WC+1:length(local_WC.res.x_pred)*(data.WC+data.L)),input_A.P(length(local_WC.res.x_pred)*(data.WC+data.L)+1:length(local_WC.res.x_pred)));

%%Zeitmatrix
% for k = 1:length(input_WC.v)-n*states-n-1
for k = 1:length(input_WC.v)-n*states-N_pred
        data.vt(k,:)=k+n-1+[0:N_pred-n];
        data.Pt(k,:)=k+n-1+[0:N_pred-n];
end
%%%Vorbereitung des Datentriplets für den Kalmanfilter
%%Datentriplet
    data.vkal = cat(1,local_WC.res.x_pred(1:end*data.WC,:),local_L.res.x_pred(end*data.WC+1:end*(data.L+data.WC),:),local_A.res.x_pred(end*(data.L+data.WC)+1:end,:));
    data.Pkal = cat(1,local_WC.res.x2_pred(1:end*data.WC,:),local_L.res.x2_pred(end*data.WC+1:end*(data.L+data.WC),:),local_A.res.x2_pred(end*(data.L+data.WC)+1:end,:));
    
%%X- & Y-Vektor
N = length(local_WC.res.x_pred); %Teilbarkeit durch 3 gegeben (s. DMDc-File)
data.X = zeros(n*states,N);
data.Xv1 = zeros(n,N);
data.XP1 = zeros(n,N);
data.Xv2 = zeros(n,N);
data.XP2 = zeros(n,N);
data.Xv3 = zeros(n,N);
data.XP3 = zeros(n,N);
data.Ypos = zeros(1,N);
data.Ypos1 = zeros(1,N);
data.Ypos2 = zeros(1,N);
data.Ypos3 = zeros(1,N);
if input == 2
    data.YE = zeros(1,N);
    data.YE1 = zeros(1,N);
    data.YE2 = zeros(1,N);
    data.YE3 = zeros(1,N);
if input == 3
    data.YE = zeros(1,N);
    data.YE1 = zeros(1,N);
    data.YE2 = zeros(1,N);
    data.YE3 = zeros(1,N);
    data.Ya = zeros(1,N);
    data.Ya1 = zeros(1,N);
    data.Ya2 = zeros(1,N);
    data.Ya3 = zeros(1,N);
elseif input == 4
    data.YE = zeros(1,N);
    data.YE1 = zeros(1,N);
    data.YE2 = zeros(1,N);
    data.YE3 = zeros(1,N);
    data.Ya = zeros(1,N);
    data.Ya1 = zeros(1,N);
    data.Ya2 = zeros(1,N);
    data.Ya3 = zeros(1,N);
    data.Yadot = zeros(1,N);
    data.Yadot1 = zeros(1,N);
    data.Yadot2 = zeros(1,N);
    data.Yadot3 = zeros(1,N);
end
end

for j = 1:n
    data.Xv1(j,:) = input_WC.v(j:N-1+j,:)';
    data.Xv2(j,:) = input_L.v(j:N-1+j,:)';
    data.Xv3(j,:) = input_A.v(j:N-1+j,:)';
    
    data.XP1(j,:) = input_WC.P(j:N-1+j,:)';
    data.XP2(j,:) = input_L.P(j:N-1+j,:)';
    data.XP3(j,:) = input_A.P(j:N-1+j,:)';
    
    data.Ypos1 = input_WC.pos(1:N);
    data.Ypos2 = input_L.pos(1:N);
    data.Ypos3 = input_A.pos(1:N);
if input == 2
    data.Ypos1 = input_WC.pos(1:N)';
    data.Ypos2 = input_L.pos(1:N)';
    data.Ypos3 = input_A.pos(1:N)';
    data.YE1 = input_WC.E(1:N)';
    data.YE2 = input_L.E(1:N)';
    data.YE3 = input_A.E(1:N)';
elseif input == 3
    data.Ypos1 = input_WC.pos(1:N)';
    data.Ypos2 = input_L.pos(1:N)';
    data.Ypos3 = input_A.pos(1:N)';
    data.YE1 = input_WC.E(1:N)';
    data.YE2 = input_L.E(1:N)';
    data.YE3 = input_A.E(1:N)';
    data.Ya1 = input_WC.adot(1:N)';
    data.Ya2 = input_L.adot(1:N)';
    data.Ya3 = input_A.adot(1:N)';
elseif input == 4
    data.Ypos1 = input_WC.pos(1:N)';
    data.Ypos2 = input_L.pos(1:N)';
    data.Ypos3 = input_A.pos(1:N)';
    data.YE1 = input_WC.E(1:N)';
    data.YE2 = input_L.E(1:N)';
    data.YE3 = input_A.E(1:N)';
    data.Ya1 = input_WC.a(1:N)';
    data.Ya2 = input_L.a(1:N)';
    data.Ya3 = input_A.a(1:N)';
    data.Ya1 = input_WC.adot(1:N)';
    data.Ya2 = input_L.adot(1:N)';
    data.Ya3 = input_A.adot(1:N)';
end
end
    data.Xv = cat(2,data.Xv1(:,1:end*data.WC),data.Xv2(:,end*data.WC+1:end*(data.L+data.WC)),data.Xv3(:,end*(data.L+data.WC)+1:end));
    data.XP = cat(2,data.XP1(:,1:end*data.WC),data.XP2(:,end*data.WC+1:end*(data.L+data.WC)),data.XP3(:,end*(data.L+data.WC)+1:end));
    data.Ypos = cat(2,data.Ypos1(:,1:end*data.WC),data.Ypos2(:,end*data.WC+1:end*(data.L+data.WC)),data.Ypos3(:,end*(data.L+data.WC)+1:end));
    if input == 2
        data.Ypos = cat(2,data.Ypos1(:,1:end*data.WC),data.Ypos2(:,end*data.WC+1:end*(data.L+data.WC)),data.Ypos3(:,end*(data.L+data.WC)+1:end));
        data.YE = cat(2,data.YE1(:,1:end*data.WC),data.YE2(:,end*data.WC+1:end*(data.L+data.WC)),data.YE3(:,end*(data.L+data.WC)+1:end));
    elseif input == 3
        data.Ypos = cat(2,data.Ypos1(:,1:end*data.WC),data.Ypos2(:,end*data.WC+1:end*(data.L+data.WC)),data.Ypos3(:,end*(data.L+data.WC)+1:end));
        data.YE = cat(2,data.YE1(:,1:end*data.WC),data.YE2(:,end*data.WC+1:end*(data.L+data.WC)),data.YE3(:,end*(data.L+data.WC)+1:end));
        data.Ya = cat(2,data.Ya1(:,1:end*data.WC),data.Ya2(:,end*data.WC+1:end*(data.L+data.WC)),data.Ya3(:,end*(data.L+data.WC)+1:end));
    elseif input == 4
        data.Ypos = cat(2,data.Ypos1(:,1:end*data.WC),data.Ypos2(:,end*data.WC+1:end*(data.L+data.WC)),data.Ypos3(:,end*(data.L+data.WC)+1:end));
        data.YE = cat(2,data.YE1(:,1:end*data.WC),data.YE2(:,end*data.WC+1:end*(data.L+data.WC)),data.YE3(:,end*(data.L+data.WC)+1:end));
        data.Ya = cat(2,data.Ya1(:,1:end*data.WC),data.Ya2(:,end*data.WC+1:end*(data.L+data.WC)),data.Ya3(:,end*(data.L+data.WC)+1:end));
        data.Yadot = cat(2,data.Yadot1(:,1:end*data.WC),data.Yadot2(:,end*data.WC+1:end*(data.L+data.WC)),data.Yadot3(:,end*(data.L+data.WC)+1:end));
    end
for ii = 1:n
    data.X(states*(ii-1)+1,:) = data.Xv(ii,:);
    data.X(states*(ii-1)+2,:) = data.XP(ii,:);
end

end