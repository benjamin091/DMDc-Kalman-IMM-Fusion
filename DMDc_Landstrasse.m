function res_L = DMDc_Landstrasse(n,N_pred,input,input_L)
    N = length(input_L.v);                  % Simulation time
    state = 2;
% Load and prepare data 
% ------------------------------------------------------------------------
%%%Dimensionierung von X und Y mit/ohne Regressoren
    X_L = zeros(n*state,N-n+1);
    y = zeros(state,length(X_L));

%%%Befüllen der Hilfsmatrizen y und z mit/ohne
%%%Regressoren

    if input == 4
        y = [round(input_L.v,4),round(input_L.P,4)]';
        z = [round(input_L.pos(1:end-n+1),4),round(input_L.E(1:end-n+1),4),round(input_L.a(1:end-n+1),4),round(input_L.adot(1:end-n+1),4)]';
    elseif input == 3
        y = [round(input_L.v,4),round(input_L.P,4)]';
        z = [round(input_L.pos(1:end-n+1),4),round(input_L.E(1:end-n+1),4),round(input_L.a(1:end-n+1),4)]';
    elseif input == 2
        y = [round(input_L.v,4),round(input_L.P,4)]';
        z = [round(input_L.pos(1:end-n+1),4),round(input_L.E(1:end-n+1),4)]';
    elseif input == 1
        y = [round(input_L.v,4),round(input_L.P,4)]';
        z = [round(input_L.pos(1:end-n+1),4)]';
    end
%%%Übertragen der Einträge aus den Hilfsmatrizen
    for ii = 1:n
        X_L((ii-1)*state+1:ii*state,:) = y(:,ii:(N-n+ii));
    end
    Y_L = z; 

%%%DMDc
    X1N1 = X_L(:,1:end-1); %X
    X2N = X_L(:,2:end); %X'
    Y1N1 = Y_L(:,1:end-1); %Y
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
    AO_L = X2N*VO/SO*UO1';
    BO_L = X2N*VO/SO*UO2';
    
    %SVD for the output space
    %=> U_dacherl, Sigma_dacherl, V_dacherl
    [UX,SX,VX] = svd(X2N,'econ');
%     [UX,SX,VX] = svd(X2N);
    
    %determination of A_tilde and B_tilde
    %A_tilde = U_dacherl^* *X'*V_tilde*Sigma_tilde^-1*U1_tilde^* *U_dacherl
    %= U_dacherl^* *A_line*U_dacherl
    %B_tilde = U_dacherl^* *X'*V_tilde*Sigma_tilde^-1*U2_tilde^*
    %= U_dacherl^* *B_line
    A_L = UX'*X2N*VO/SO*UO1'*UX;
    B_L = UX'*X2N*VO/SO*UO2';
    
    x_pred = zeros(n*state,N_pred);
    u_pred = zeros(input,N_pred);

    res_L=struct;
    res_L.x_pred=zeros(N-n*state-N_pred,N_pred);
    res_L.x2_pred=zeros(N-n*state-N_pred,N_pred);
       
    for ii=1:1:N-n*state-N_pred
        %wird an die erste Spalte von x_pred
        %es werden N_pred Einträge von Y an u_pred übergeben und um eins weiter nach rechts versetzt 
        x_pred(:,1) = X_L(:,ii);
        u_pred(:,1:N_pred) = Y_L(:,ii:(N_pred+ii-1));
        %über die Systemmatrizen A und B wird die die kk+1-te Spalte 
        %der Systemzustände prädiziert
        for kk=1:N_pred-1
%               x_pred(:,kk+1) = A*x_pred(:,kk)+B*u_pred(:,kk);
              x_pred(:,kk+1) = AO_L*x_pred(:,kk)+BO_L*u_pred(:,kk);      
        end
        
        %die prädizierten Zustände werden an die Größe res.x_pred übergeben
        res_L.x_pred(ii,:) = x_pred(1,:);
        res_L.t_pred(ii,:) = ii+n-1+[0:N_pred];
        res_L.x2_pred(ii,:) = x_pred(2,:);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Teilbarkeit durch 3 gewährleisten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while mod(length(res_L.x_pred),3) ~= 0
        res_L.x_pred = res_L.x_pred(1:end-1,:);
        res_L.x2_pred = res_L.x2_pred(1:end-1,:);
    end
    res_L.A = AO_L;
    res_L.B = BO_L;
    res_L.X = X_L;
    res_L.Y = Y_L;
end