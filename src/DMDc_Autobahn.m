function res_A = DMDc_Autobahn(n,N_pred,input,input_A)
    N = length(input_A.v);                  % Simulation time
    state = 2;
% Load and prepare data 
% ------------------------------------------------------------------------
%%%Dimensionierung von X und Y mit/ohne Regressoren
    X_A = zeros(n*state,N-n+1);
    y = zeros(state,length(X_A));

%%%Befüllen der Hilfsmatrizen y und z mit/ohne
%%%Regressoren
    if input == 4
        y = [round(input_A.v,4),round(input_A.P,4)]';
        z = [round(input_A.pos(1:end-n+1),4),round(input_A.E(1:end-n+1),4),round(input_A.a(1:end-n+1),4),round(input_A.adot(1:end-n+1),4)]';
    elseif input == 3
        y = [round(input_A.v,4),round(input_A.P,4)]';
        z = [round(input_A.pos(1:end-n+1),4),round(input_A.E(1:end-n+1),4),round(input_A.a(1:end-n+1),4)]';
    elseif input == 2
        y = [round(input_A.v,4),round(input_A.P,4)]';
        z = [round(input_A.pos(1:end-n+1),4),round(input_A.E(1:end-n+1),4)]';
    elseif input == 1
        y = [round(input_A.v,4),round(input_A.P,4)]';
        z = [round(input_A.pos(1:end-n+1),4)]';
    end

%%%Übertragen der Einträge aus den Hilfsmatrizen
    for ii = 1:n
        X_A((ii-1)*state+1:ii*state,:) = y(:,ii:(N-n+ii));
    end
    Y_A = z; 

%%%DMDc
    X1N1 = X_A(:,1:end-1); %X
    X2N = X_A(:,2:end); %X'
    Y1N1 = Y_A(:,1:end-1); %Y
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
    AO_A = X2N*VO/SO*UO1';
    BO_A = X2N*VO/SO*UO2';
    
    %SVD for the output space
    %=> U_dacherl, Sigma_dacherl, V_dacherl
    [UX,SX,VX] = svd(X2N,'econ');
%     [UX,SX,VX] = svd(X2N);
    
    %determination of A_tilde and B_tilde
    %A_tilde = U_dacherl^* *X'*V_tilde*Sigma_tilde^-1*U1_tilde^* *U_dacherl
    %= U_dacherl^* *A_line*U_dacherl
    %B_tilde = U_dacherl^* *X'*V_tilde*Sigma_tilde^-1*U2_tilde^*
    %= U_dacherl^* *B_line
    A_A = UX'*X2N*VO/SO*UO1'*UX;
    B_A = UX'*X2N*VO/SO*UO2';
    
    x_pred = zeros(n*state,N_pred);
    u_pred = zeros(input,N_pred);

    res_A=struct;
    res_A.x_pred=zeros(N-n*state-N_pred,N_pred);
    res_A.x2_pred=zeros(N-n*state-N_pred,N_pred);
       
    for ii=1:1:N-n*state-N_pred
        %wird an die erste Spalte von x_pred
        %es werden N_pred Einträge von Y an u_pred übergeben und um eins weiter nach rechts versetzt 
        x_pred(:,1) = X_A(:,ii);
        u_pred(:,1:N_pred) = Y_A(:,ii:(N_pred+ii-1));
        %über die Systemmatrizen A und B wird die die kk+1-te Spalte 
        %der Systemzustände prädiziert
        for kk=1:N_pred-1
%               x_pred(:,kk+1) = A_A*x_pred(:,kk)+B_A*u_pred(:,kk);
              x_pred(:,kk+1) = AO_A*x_pred(:,kk)+BO_A*u_pred(:,kk);      
        end
        
        %die prädizierten Zustände werden an die Größe res.x_pred übergeben
        res_A.x_pred(ii,:) = x_pred(1,:);
        res_A.t_pred(ii,:) = ii+n-1+[0:N_pred];
        res_A.x2_pred(ii,:) = x_pred(2,:);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Teilbarkeit durch 3 gewährleisten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     while mod(length(res_A.x_pred),3) ~= 0
        res_A.x_pred = res_A.x_pred(1:end-1,:);
        res_A.x2_pred = res_A.x2_pred(1:end-1,:);  
    end
    res_A.A = AO_A;
    res_A.B = BO_A;
    res_A.X = X_A;
    res_A.Y = Y_A;
end