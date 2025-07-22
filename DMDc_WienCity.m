function res_WC = DMDc_WienCity(n,N_pred,input,input_WC)
    N = length(input_WC.v);                  % Simulation time
    state = 2;
% Load and prepare data 
% ------------------------------------------------------------------------
% Dimensionierung von X und Y 
X_WC = zeros(n*state,N-n+1);
y = zeros(state,length(X_WC));

%%%Befüllen der Hilfsmatrizen y und z 
if input == 4
    y = [round(input_WC.v,4),round(input_WC.P,4)]';
    z = [round(input_WC.pos(1:end-n+1),4),round(input_WC.E(1:end-n+1),4),round(input_WC.a(1:end-n+1),4),round(input_WC.adot(1:end-n+1),4)]';
elseif input == 3
    y = [round(input_WC.v,4),round(input_WC.P,4)]';
    z = [round(input_WC.pos(1:end-n+1),4),round(input_WC.E(1:end-n+1),4),round(input_WC.a(1:end-n+1),4)]';
elseif input == 2
    y = [round(input_WC.v,4),round(input_WC.P,4)]';
    z = [round(input_WC.pos(1:end-n+1),4),round(input_WC.E(1:end-n+1),4)]';
elseif input == 1
    y = [round(input_WC.v,4),round(input_WC.P,4)]';
    z = [round(input_WC.pos(1:end-n+1),4)]';
end

%Übertragen der Einträge aus den Hilfsmatrizen
for ii = 1:n
    X_WC((ii-1)*state+1:ii*state,:) = y(:,ii:(N-n+ii));
end
Y_WC = z; 

%%%DMDc
    X1N1 = X_WC(:,1:end-1); %X
    X2N = X_WC(:,2:end); %X'
    Y1N1 = Y_WC(:,1:end-1); %Y
    O = [X1N1;Y1N1]; %Omega_Groß
    
    %SVD on the augmented data; input space
    %=> U_tilde,Sigma_tilde,V_tilde
    [UO,SO,VO] = svd(O,'econ'); 
%     [UO,SO,VO] = svd(O); 

    %decomposition of U_tilde into U1_tilde and U2_tilde
    %dim(U1_tilde) = n*states x p (=truncation value)
    %dim(U2_tilde) = input x p

UO1 = UO(1:state*n,1:end); 
UO2 = UO((state*n)+1:end,1:end); 
    %determination of A_line and B_line with the decomposed Ui_tilde
    %[A_line, B_line] = [X'*V_tilde*Sigma_tilde^-1*U1_tilde^*,
    %X'*V_tilde*Sigma_tilde^-1*U2_tilde*]
AO_WC = X2N*VO/SO*UO1';
BO_WC = X2N*VO/SO*UO2';
    
    %SVD for the output space
    %=> U_dacherl, Sigma_dacherl, V_dacherl
[UX,SX,VX] = svd(X2N,'econ');
%     [UX,SX,VX] = svd(X2N);
    
    %determination of A_tilde and B_tilde
    %A_tilde = U_dacherl^* *X'*V_tilde*Sigma_tilde^-1*U1_tilde^* *U_dacherl
    %= U_dacherl^* *A_line*U_dacherl
    %B_tilde = U_dacherl^* *X'*V_tilde*Sigma_tilde^-1*U2_tilde^*
    %= U_dacherl^* *B_line
%     A = UX'*AO*UX;
%     B = UX'*BO;
A_WC = UX'*X2N*VO/SO*UO1'*UX;
B_WC = UX'*X2N*VO/SO*UO2';
    
x_pred = zeros(n*state,N_pred);
u_pred = zeros(input,N_pred);

res_WC=struct;
res_WC.x_pred=zeros(N-n*state-N_pred,N_pred);
res_WC.x2_pred=zeros(N-n*state-N_pred,N_pred);

for ii=1:1:N-n*state-N_pred
    x_pred(:,1) = X_WC(:,ii);
    u_pred(:,1:N_pred) = Y_WC(:,ii:(N_pred+ii-1));
    for kk=1:N_pred-1
%               x_pred(:,kk+1) = A_WC*x_pred(:,kk)+B_WC*u_pred(:,kk);
        x_pred(:,kk+1) = AO_WC*x_pred(:,kk)+BO_WC*u_pred(:,kk);      
    end
    res_WC.x_pred(ii,:) = x_pred(1,:);
    res_WC.t_pred(ii,:) = ii+n-1+[0:N_pred];
    res_WC.x2_pred(ii,:) = x_pred(2,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Teilbarkeit durch 3 gewährleisten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while mod(length(res_WC.x_pred),3) ~= 0
        res_WC.x_pred = res_WC.x_pred(1:end-1,:);
        res_WC.x2_pred = res_WC.x2_pred(1:end-1,:);
    end
    res_WC.A = AO_WC;
    res_WC.B = BO_WC;
    res_WC.X = X_WC;
    res_WC.Y = Y_WC;
end