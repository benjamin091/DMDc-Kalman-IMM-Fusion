%% DMDc
% 
% ########################################################################
% Beschreibung:         correlation is a measure of similarity between two signals, and convolution is a measure of effect of one signal on the other.
% Version:              V1
% erstellt am:          27.09.2021
% letzte Änderung am:   28.29.2021
% Änderungen:           -)  Erstversion (24.03.2021)
close all; 
clear all;
clc           
% Anmerkungen:          -) ...
% Quellen:              
% ########################################################################

% Cleaning 
% ------------------------------------------------------------------------

%%Globale Variabeln
t_abtast = 1;
n = 8;
N_pred = 20;
input = 2;
states = 2;

const = struct;
const.m = 1746;
const.g = 9.81;
const.cr = 0.01;
const.Av = 8;
const.cx = 0.35;
const.rho = 1.2;
const.fall = 3;

%Anteile an Strassentypen - insgesamt 3 Anteile
data_WC = 1;
data_L = 1;
data_A = 3 - data_WC - data_L;

train_val = true; %%%trainingsdaten true/validierungsdaten false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Daten aus files laden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%WienCity
input_WC = struct;
if train_val == true
    %%Trainingsdaten - VW Passat SDriver
    load('Passat_SDriver_WienCity.mat')
    input_WC.v = Passat_SDriver_WienCity.velocity_kmh;
    input_WC.a = Passat_SDriver_WienCity.acc;
    input_WC.t = Passat_SDriver_WienCity.t;
    input_WC.E = Passat_SDriver_WienCity.energy;
    input_WC.slope = Passat_SDriver_WienCity.slope;
    input_WC.lat = Passat_SDriver_WienCity.lat;
    input_WC.lon = Passat_SDriver_WienCity.lon;
else
%%Validierungsdaten - VW FC SDriver
    load('Passat_FCDriver_WienCity.mat')
    input_WC.v = Passat_FCDriver_WienCity.velocity_kmh;
    input_WC.a = Passat_FCDriver_WienCity.acc;
    input_WC.t = Passat_FCDriver_WienCity.t;
    input_WC.E = Passat_FCDriver_WienCity.energy;
    input_WC.slope = Passat_FCDriver_WienCity.slope;
    input_WC.lat = Passat_FCDriver_WienCity.lat;
    input_WC.lon = Passat_FCDriver_WienCity.lon;
end

    input_WC.v = input_WC.v(~isnan(input_WC.v));
    input_WC.a = input_WC.a(~isnan(input_WC.a));
    input_WC.E = input_WC.E(~isnan(input_WC.E));
    input_WC.slope = input_WC.slope(~isnan(input_WC.slope));
    input_WC.lat = input_WC.lat(~isnan(input_WC.lat));
    input_WC.lon = input_WC.lon(~isnan(input_WC.lon));
    
    input_WC.v = smooth(input_WC.v);
    input_WC.a = smooth(input_WC.a);
    input_WC.E = smooth(input_WC.E);
    input_WC.slope = smooth(input_WC.slope);
    input_WC.lat = smooth(input_WC.lat);
    input_WC.lon = smooth(input_WC.lon);

%%%Landstrasse
input_L = struct;
if train_val == true
    %%Trainingsdaten - VW Passat SDriver
    load('Passat_SDriver_Landstrasse.mat');
    input_L.v = Passat_SDriver_Landstrasse.velocity_kmh;
    input_L.a = Passat_SDriver_Landstrasse.acc;
    input_L.t = Passat_SDriver_Landstrasse.t;
    input_L.E = Passat_SDriver_Landstrasse.energy;
    input_L.slope = Passat_SDriver_Landstrasse.slope;
    input_L.lat = Passat_SDriver_Landstrasse.lat;
    input_L.lon = Passat_SDriver_Landstrasse.lon;
else
    %%Validierungsdaten - VW FC SDriver
    load('Passat_FCDriver_Landstrasse.mat')
    input_L.v = Passat_FCDriver_Landstrasse.velocity_kmh;
    input_L.a = Passat_FCDriver_Landstrasse.acc;
    input_L.t = Passat_FCDriver_Landstrasse.t;
    input_L.E = Passat_FCDriver_Landstrasse.energy;
end
    input_L.v = input_L.v(~isnan(input_L.v));
    input_L.a = input_L.a(~isnan(input_L.a));
    input_L.E = input_L.E(~isnan(input_L.E));
    input_L.slope = input_L.slope(~isnan(input_L.slope));
    input_L.lat = input_L.slope(~isnan(input_L.lat));
    input_L.lon = input_L.lon(~isnan(input_L.lon));
    
    input_L.v = smooth(input_L.v);
    input_L.a = smooth(input_L.a);
    input_L.E = smooth(input_L.E);
    input_L.slope = smooth(input_L.slope);
    input_L.lat = smooth(input_L.lat);
    input_L.lon = smooth(input_L.lon);
  
    
%%%Autobahn
input_A = struct;
if train_val == true
    %%Trainingsdaten - VW Passat SDriver
    load('Passat_SDriver_Autobahntraffic.mat')
    input_A.v = Passat_SDriver_Autobahntraffic.velocity_kmh;
    input_A.a = Passat_SDriver_Autobahntraffic.acc;
    input_A.t = Passat_SDriver_Autobahntraffic.t;
    input_A.E = Passat_SDriver_Autobahntraffic.energy;
    input_A.slope = Passat_SDriver_Autobahntraffic.slope;
    input_A.lat = Passat_SDriver_Autobahntraffic.lat;
    input_A.lon = Passat_SDriver_Autobahntraffic.lon;
else
%%Validierungsdaten - VW FC SDriver
    load('Passat_FCDriver_Autobahn.mat')
    input_A.v = Passat_FCDriver_Autobahn.velocity_kmh;
    input_A.a = Passat_FCDriver_Autobahn.acc;
    input_A.t = Passat_FCDriver_Autobahn.t;
    input_A.E = Passat_FCDriver_Autobahn.energy;
end

    input_A.v = input_A.v(~isnan(input_A.v));
    input_A.a = input_A.a(~isnan(input_A.a));
    input_A.E = input_A.E(~isnan(input_A.E));
    input_A.slope = input_A.slope(~isnan(input_A.slope));
    input_A.lat = input_A.lat(~isnan(input_A.lat));
    input_A.lon = input_A.lon(~isnan(input_A.lon));
    
    input_A.v = smooth(input_A.v);
    input_A.a = smooth(input_A.a);
    input_A.E = smooth(input_A.E);
    input_A.slope = smooth(input_A.slope);
    input_A.lat = smooth(input_A.lat);
    input_A.lon = smooth(input_A.lon);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Daten samplen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_data = min([length(input_WC.v),length(input_L.v),length(input_A.v)]);


    input_WC.v = input_WC.v(1:t_abtast/.5:min_data);
    input_WC.a = input_WC.a(1:t_abtast/.5:min_data);
    input_WC.t = input_WC.t(1:t_abtast/.5:min_data);
    input_WC.E = input_WC.E(1:t_abtast/.5:min_data);
    input_WC.slope = input_WC.slope(1:t_abtast/.5:min_data);
    input_WC.lat = input_WC.lat(1:t_abtast/.5:min_data);
    input_WC.lon = input_WC.lon(1:t_abtast/.5:min_data);

    input_L.v = input_L.v(1:t_abtast/.5:min_data);
    input_L.a = input_L.a(1:t_abtast/.5:min_data);
    input_L.t = input_L.t(1:t_abtast/.5:min_data);
    input_L.E = input_L.E(1:t_abtast/.5:min_data);
    input_L.slope = input_L.slope(1:t_abtast/.5:min_data);
    input_L.lat = input_L.lat(1:t_abtast/.5:min_data);
    input_L.lon = input_L.lon(1:t_abtast/.5:min_data);
   
    
    input_A.v = input_A.v(1:t_abtast/.5:min_data);
    input_A.a = input_A.a(1:t_abtast/.5:min_data);
    input_A.t = input_A.t(1:t_abtast/.5:min_data);
    input_A.E = input_A.E(1:t_abtast/.5:min_data);
    input_A.slope = input_A.slope(1:t_abtast/.5:min_data);
    input_A.lat = input_A.lat(1:t_abtast/.5:min_data);
    input_A.lon = input_A.lon(1:t_abtast/.5:min_data);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%abgeleitete Parameter bestimmen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wien City
input_WC.t_delta = input_WC.t(2)-input_WC.t(1);
if const.fall == 1
    input_WC.P = [0; diff(input_WC.E)/input_WC.t_delta]; 
elseif const.fall == 2
    input_WC.Froll = const.m*const.g*const.cr*cos(deg2rad(input_WC.slope));
    input_WC.Fslope = const.m*const.g*sin(deg2rad(input_WC.slope));
    input_WC.Fdrag = 1/2*const.Av*const.cx*const.rho*(input_WC.v/3.6).^2;
    input_WC.Fres = input_WC.Froll + input_WC.Fslope + input_WC.Fdrag;
    input_WC.P = (const.m*input_WC.a + input_WC.Fres).*input_WC.v/3.6';
elseif const.fall == 3
    input_WC.slope = 0;
    input_WC.Froll = const.m*const.g*const.cr*cos(deg2rad(input_WC.slope));
    input_WC.Fslope = const.m*const.g*sin(deg2rad(input_WC.slope));
    input_WC.Fdrag = 1/2*const.Av*const.cx*const.rho*(input_WC.v/3.6).^2;
    input_WC.Fres = input_WC.Froll + input_WC.Fslope + input_WC.Fdrag;
    input_WC.P = (const.m*input_WC.a + input_WC.Fres).*input_WC.v/3.6';
end             
input_WC.adot = [0; diff(input_WC.a)/input_WC.t_delta];
input_WC.pos = sqrt(input_WC.lat.^2+input_WC.lon.^2);
input_WC.E = [0; diff(input_WC.E)];

%%Landstrasse
input_L.t_delta = input_L.t(2)-input_L.t(1);
if const.fall == 1
    input_L.P = [0; diff(input_L.E)/input_L.t_delta]; 
elseif const.fall ==2
    input_L.Froll = const.m*const.g*const.cr*cos(deg2rad(input_L.slope));
    input_L.Fslope = const.m*const.g*sin(deg2rad(input_L.slope));
    input_L.Fdrag = 1/2*const.Av*const.cx*const.rho*(input_L.v/3.6).^2;
    input_L.Fres = input_L.Froll + input_L.Fslope + input_L.Fdrag;
    input_L.P = (const.m*input_L.a + input_L.Fres).*input_L.v/3.6';
elseif const.fall == 3
    input_L.slope = 0;
    input_L.Froll = const.m*const.g*const.cr*cos(deg2rad(input_L.slope));
    input_L.Fslope = const.m*const.g*sin(deg2rad(input_L.slope));
    input_L.Fdrag = 1/2*const.Av*const.cx*const.rho*(input_L.v/3.6).^2;
    input_L.Fres = input_L.Froll + input_L.Fslope + input_L.Fdrag;
    input_L.P = (const.m*input_L.a + input_L.Fres).*input_L.v/3.6';
end
input_L.adot = [0; diff(input_L.a)/input_L.t_delta];
input_L.pos = sqrt(input_L.lat.^2+input_L.lon.^2);
input_L.E = [0; diff(input_L.E)];

%%Autobahn
input_A.t_delta = input_A.t(2)-input_A.t(1);
if const.fall == 1
    input_A.P = [0; diff(input_A.E)/input_A.t_delta]; 
elseif const.fall == 2
    input_A.Froll = const.m*const.g*const.cr*cos(deg2rad(input_A.slope));
    input_A.Fslope = const.m*const.g*sin(deg2rad(input_A.slope));
    input_A.Fdrag = 1/2*const.Av*const.cx*const.rho*(input_A.v/3.6).^2;
    input_A.Fres = input_A.Froll + input_A.Fslope + input_A.Fdrag;
    input_A.P = (const.m*input_A.a + input_A.Fres).*input_A.v/3.6';
elseif const.fall == 3
    input_A.slope = 0;
    input_A.Froll = const.m*const.g*const.cr*cos(deg2rad(input_A.slope));
    input_A.Fslope = const.m*const.g*sin(deg2rad(input_A.slope));
    input_A.Fdrag = 1/2*const.Av*const.cx*const.rho*(input_A.v/3.6).^2;
    input_A.Fres = input_A.Froll + input_A.Fslope + input_A.Fdrag;
    input_A.P = (const.m*input_A.a + input_A.Fres).*input_A.v/3.6';

end              
input_A.adot = [0; diff(input_A.a)/input_A.t_delta];  
input_A.pos = sqrt(input_A.lat.^2+input_A.lon.^2);
input_A.E = [0; diff(input_A.E)];

%%%%%%%%%%%%%%%%%%%%%%%%
%%DMDc ausführen
%%%%%%%%%%%%%%%%%%%%%%%
res_WC = DMDc_WienCity(n,N_pred,input,input_WC);
res_L = DMDc_Landstrasse(n,N_pred,input,input_L);
res_A = DMDc_Autobahn(n,N_pred,input,input_A); 

%%%%%%%%%%%%%%%%%%%%%%%%
%%Stability Plot
%%%%%%%%%%%%%%%%%%%%%%%%
% figure,plot(eig(res_WC.A),'bo');zgrid
% legend('A-Matrix Wien City','Interpreter','latex','FontSize',10);
% axis equal

% figure,plot(eig(res_L.A),'bo');zgrid
% legend6 = legend('A-Matrix Landstrasse','Interpreter','latex','FontSize',10);
% axis equal

% figure,plot(eig(res_A.A),'bo');zgrid
% legend5 = legend('A-Matrix Autobahn','Interpreter','latex','FontSize',10);
% axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%
%%Kalman-Filter
%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wien City
local_WC = struct;
% local_WC.Q = diag([1 .1]);
% local_WC.R = diag([.1 100]);
local_WC.Q = diag([1]);
local_WC.R = diag([1 1]);
local_WC.N = 0;
local_WC.T = t_abtast;
local_WC.res = res_WC;

kal_WC = Kalman_DMDc_WienCity(local_WC,n,N_pred,input);

%%Landstrasse
local_L = struct;
% local_L.Q = diag([1 .1]);
% local_L.R = diag([.1 100]);
local_L.Q = diag([1]);
local_L.R = diag([1 1]);
local_L.N = 0;
local_L.T = t_abtast;
local_L.res = res_L;

kal_L = Kalman_DMDc_Landstrasse(local_L,n,N_pred,input);

%%Autobahn
local_A = struct;
% local_A.Q = diag([1 .1]); 
% local_A.R = diag([.1 10]);
local_A.Q = diag([1]); 
local_A.R = diag([1 1]);
local_A.N = 0;
local_A.T = t_abtast;
local_A.res = res_A;

kal_A = Kalman_DMDc_Autobahn(local_A,n,N_pred,input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Daten für IMM und Referenzkalmanfilter vorbereiten 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = data_fct_DMDc(input_WC,input_L,input_A,local_WC,local_L,local_A,N_pred,n,states,input,data_WC,data_L,data_A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%IMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markov = struct;
markov.p11 = .95;
markov.p22 = .95;
markov.p33 = .95;
markov.m = 1;
mu_hat = [.33; .33; .33];
[IMM,P_IMM] = IMM_DMDc(n,N_pred,markov,mu_hat,local_WC,local_L,local_A,kal_WC,kal_L,kal_A,data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Referenzkalmanfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_WC_test = Kalman_DMDc_WienCity_Test(local_WC,n,N_pred,input,data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Keine Datenmischung - Kalman vs IMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure('units','normalized','outerposition',[0.0 0.0 1.0 1.0])
%     ax1 = subplot(2,1,1);
%     plot(data.vref,"k",'LineWidth',1.2)
%     hold on
%     for pl=1:15:length(data.vref)
%         if mod(pl,2)
%             plot(data.vt(pl,:),kal_WC.sim.vplot(pl,:),"r:",'LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.vt(pl,1), kal_WC.sim.vplot(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.vt(pl,end), kal_WC.sim.vplot(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%         else
%             plot(data.vt(pl,:),kal_WC.sim.vplot(pl,:),'g:','LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.vt(pl,1), kal_WC.sim.vplot(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.vt(pl,end), kal_WC.sim.vplot(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
% 
%         end
%     end
%     legend2 = legend('v in km/h','v in km/h estimated','FontSize',10);
%     title('Kalman (current): 8 Regression steps, 12 Prediction time steps, WienCity data');
%     set(legend2,'Orientation','horizontal','FontSize',10);
%     grid on
% 
%     ax2 = subplot(2,1,2);
%     plot(data.vref,"k",'LineWidth',1.2)
%     hold on
%     for pl=1:15:length(data.vref)
%         if mod(pl,2)
%             plot(data.vt(pl,:),IMM.v(pl,:),"r:",'LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.vt(pl,1), IMM.v(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.vt(pl,end), IMM.v(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%         else
%             plot(data.vt(pl,:),IMM.v(pl,:),'g:','LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.vt(pl,1), IMM.v(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.vt(pl,end), IMM.v(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
% 
%         end
%     end
%     legend2 = legend('v in km/h','v in km/h estimated','FontSize',10);
%     title('IMM: 8 Regression steps, 12 Prediction time steps, WienCity data');
%     set(legend2,'Orientation','horizontal','FontSize',10);
%     grid on
%     xlabel('time in s');
%     linkaxes([ax1,ax2],'x')
%     sgtitle('Kalman vs IMM - WienCity');
% 
%     figure('units','normalized','outerposition',[0.0 0.0 1.0 1.0])
%     ax1 = subplot(2,1,1);
%     plot(data.Pref,"k",'LineWidth',1.2)
%     hold on
%     for pl=1:15:length(data.Pref)
%         if mod(pl,2)
%             plot(data.Pt(pl,:),kal_WC.sim.Pplot(pl,:),"r:",'LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.Pt(pl,1), kal_WC.sim.Pplot(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.Pt(pl,end), kal_WC.sim.Pplot(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%         else
%             plot(data.Pt(pl,:),kal_WC.sim.Pplot(pl,:),'g:','LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.Pt(pl,1), kal_WC.sim.Pplot(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.Pt(pl,end), kal_WC.sim.Pplot(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
% 
%         end
%     end
%     grid on
%     legend2 = legend('P in Watt','P in Watt estimated','FontSize',10);
%     title('Kalman (current): 8 Regression steps, 12 Prediction time steps; WienCity data');
%     set(legend2,'Orientation','horizontal','FontSize',10);
%     
%     ax2 = subplot(2,1,2);
%     plot(data.Pref,"k",'LineWidth',1.2)
%     hold on
%     for pl=1:15:length(data.Pref)
%         if mod(pl,2)
%             plot(data.Pt(pl,:),IMM.P(pl,:),"r:",'LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.Pt(pl,1), IMM.P(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.Pt(pl,end), IMM.P(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%         else
%             plot(data.Pt(pl,:),IMM.P(pl,:),'g:','LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.Pt(pl,1), IMM.P(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.Pt(pl,end), IMM.P(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
% 
%         end
%     end
%     grid on
%     legend2 = legend('P in Watt','P in Watt estimated','FontSize',10);
%     title('IMM: 8 Regression steps, 12 Prediction time steps, WienCity data');
%     set(legend2,'Orientation','horizontal','FontSize',10);
%     xlabel('time in s');
%     linkaxes([ax1,ax2],'x')
%     sgtitle('Kalman vs IMM - WienCity');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Datentriplet: Kalman vs IMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('units','normalized','outerposition',[0.0 0.0 1.0 1.0])
    ax1 = subplot(2,1,1);
    plot(data.vref,"k",'LineWidth',1.2)
    hold on
    for pl=1:15:length(data.vref)
        if mod(pl,2)
            plot(data.vt(pl,:),sim_WC_test.v(pl,:),"r:",'LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.vt(pl,1), sim_WC_test.v(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.vt(pl,end), sim_WC_test.v(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
        else
            plot(data.vt(pl,:),sim_WC_test.v(pl,:),'g:','LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.vt(pl,1), sim_WC_test.v(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.vt(pl,end), sim_WC_test.v(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);

        end
    end
    legend2 = legend('v in km/h','v in km/h estimated','FontSize',10);
    title('Kalman (current): 8 Regression steps, 12 Prediction time steps');
    set(legend2,'Orientation','horizontal','FontSize',10);
    grid on

    ax2 = subplot(2,1,2);
    plot(data.vref,"k",'LineWidth',1.2)
    hold on
    for pl=1:15:length(data.vref)
        if mod(pl,2)
            plot(data.vt(pl,:),IMM.v(pl,:),"r:",'LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.vt(pl,1), IMM.v(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.vt(pl,end), IMM.v(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
        else
            plot(data.vt(pl,:),IMM.v(pl,:),'g:','LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.vt(pl,1), IMM.v(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.vt(pl,end), IMM.v(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);

        end
    end
    legend2 = legend('v in km/h','v in km/h estimated','FontSize',10);
    title('IMM: 8 Regression steps, 12 Prediction time steps');
    set(legend2,'Orientation','horizontal','FontSize',10);
    grid on
    xlabel('time in s');
    linkaxes([ax1,ax2],'x')
    sgtitle('Kalman vs IMM - 1/3 of each data type')
    
    figure('units','normalized','outerposition',[0.0 0.0 1.0 1.0])
    ax1 = subplot(2,1,1);
    plot(data.Pref,"k",'LineWidth',1.2)
    hold on
    for pl=1:15:length(data.Pref)
        if mod(pl,2)
            plot(data.Pt(pl,:),sim_WC_test.P(pl,:),"r:",'LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.Pt(pl,1), sim_WC_test.P(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.Pt(pl,end), sim_WC_test.P(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
        else
            plot(data.Pt(pl,:),sim_WC_test.P(pl,:),'g:','LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.Pt(pl,1), sim_WC_test.P(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.Pt(pl,end), sim_WC_test.P(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);

        end
    end
    grid on
    legend2 = legend('P in Watt','P in Watt estimated','FontSize',10);
    title('Kalman (current): 8 Regression steps, 12 Prediction time steps');
    set(legend2,'Orientation','horizontal','FontSize',10);
    
    ax2 = subplot(2,1,2);
    plot(data.Pref,"k",'LineWidth',1.2)
    hold on
    for pl=1:15:length(data.Pref)
        if mod(pl,2)
            plot(data.Pt(pl,:),IMM.P(pl,:),"r:",'LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.Pt(pl,1), IMM.P(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.Pt(pl,end), IMM.P(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
        else
            plot(data.Pt(pl,:),IMM.P(pl,:),'g:','LineWidth',1.2)
            labels =  {sprintf('$o_{%d}$', pl) };
            labelpoints(data.Pt(pl,1), IMM.P(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
            labels =  {sprintf('$x_{%d}$', pl) };
            labelpoints(data.Pt(pl,end), IMM.P(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);

        end
    end
    grid on
    legend2 = legend('P in Watt','P in Watt estimated','FontSize',10);
    title('IMM: 8 Regression steps, 12 Prediction time steps');
    set(legend2,'Orientation','horizontal','FontSize',10);
    xlabel('time in s');
    linkaxes([ax1,ax2],'x')
    sgtitle('Kalman vs IMM - 1/3 of each data type')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Datentriplet: IMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure('units','normalized','outerposition',[0.0 0.0 1.0 1.0])
%     ax1 = subplot(2,1,1);
%     plot(data.vref,"k",'LineWidth',1.2)
%     hold on
%     for pl=1:15:length(data.vref)
%         if mod(pl,2)
%             plot(data.vt(pl,:),IMM.v(pl,:),"r:",'LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.vt(pl,1), IMM.v(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.vt(pl,end), IMM.v(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%         else
%             plot(data.vt(pl,:),IMM.v(pl,:),'g:','LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.vt(pl,1), IMM.v(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.vt(pl,end), IMM.v(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
% 
%         end
%     end
%     legend2 = legend('v in km/h','v in km/h estimated','FontSize',10);
%     title('IMM: 8 Regression steps, 12 Prediction time steps');
%     set(legend2,'Orientation','horizontal','FontSize',10);
%     grid on
% 
%     ax2 = subplot(2,1,2);
%     plot(data.Pref,"k",'LineWidth',1.2)
%     hold on
%     for pl=1:15:length(data.Pref)
%         if mod(pl,2)
%             plot(data.Pt(pl,:),IMM.P(pl,:),"r:",'LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.Pt(pl,1), IMM.P(pl,1),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.Pt(pl,end), IMM.P(pl,end),labels,'interpreter','latex', 'Color', 'r','FontSize', 11);
%         else
%             plot(data.Pt(pl,:),IMM.P(pl,:),'g:','LineWidth',1.2)
%             labels =  {sprintf('$o_{%d}$', pl) };
%             labelpoints(data.Pt(pl,1), IMM.P(pl,1),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
%             labels =  {sprintf('$x_{%d}$', pl) };
%             labelpoints(data.Pt(pl,end), IMM.P(pl,end),labels,'interpreter','latex', 'Color', 'g','FontSize', 11);
% 
%         end
%     end
%     grid on
%     legend2 = legend('P in Watt','P in Watt estimated','FontSize',10);
%     title('IMM: 8 Regression steps, 12 Prediction time steps');
%     set(legend2,'Orientation','horizontal','FontSize',10);
%     xlabel('time in s');
%     linkaxes([ax1,ax2],'x')
%     sgtitle('IMM: 1/3 of each street type data');
%     