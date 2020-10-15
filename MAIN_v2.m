%% ==================================================================== %%
% Creator: Joseph N. Squeo
% Contact: joseph.squeo@uconn.edu
% Version History: v1 - 4/6/2020
% Affiliation: Univeristy of Connecticut

% Example taken from:
%[1]        S. Gillijns, O. B. Mendoza, J. Chandrasekar, B. L. R. de Moor,
% D. S. Bernstein, and A. Ridley, What is the ensemble Kalman filter and
% how well does it work?, 2006.

% Goal:
%   Improve prediciton of state variables (sootFv, temperature, ...) by
%   assimilating measurements and their uncertainty into a model, in this
%   case the two-equation soot model based on that of Lindstedt. et al.

% Function:
%   _Main_StanfordEnKF.m_ is the main MATLAB script used to initiate
%   the EnKF parameters and OpenFOAM simulation and call the EnKF script
%   called _santoro_EnKF_Stanford.m_. This script is coded to run the Santoro
%   laminar diffusion ethylene/air flame in OpenFOAM on the local Yosemite
%   Linux machine. The solver name must be specified, along with the path to
%   the case folder. Measurements for assimilation are stored into a matrix
%   in this script, where the matrix index corresponds to the cell location
%   in the mesh from OpenFOAM. This is found by using the cell centers Cx, Cy
%   and Cz in OpenFOAM which corresponds to the measurement location at some
%   radial location and height above the burner (HAB). Separate functions are
%   used to postProcess the data in OpenFOAM with the sampleDict utility for
%   plotting results at the end.

% CASE SET-UP:
%   1. Ensure the paths in this script are set to proper locations
%   2. Set up system/decomposeParDict file in the case folder
%   3. Create a file named 'callSolverParallel' which calles the solver in
%      parallel with the specified number of processors 
%        i.e.) nohup mpirun -np 20 fvJacob -parallel
%   4. Ensure start time folder is in the case folder
%   5. generatePlots.m requires there to be the baseline data as a
%   file i.e. nominal/postProcessing/sampleDict/1.93. Set the path using
%   EnkF.baseLine
%  ==================================================================== %

clear
close all
clc
addpath('~/dataAssimilation/1DheatConduction');                                    % specify path to MATLAB functions

%% === 1.Initialization ===

% Material Properties (aluminum) 
kTherm = 80;                                                               %[W / m*K] thermal conductivity
Cp = 460;                                                                  %[J / kg*K] specific heat @ constant pressure
rho = 7870;                                                                %[kg / m^3] density
rhoCp = rho .* Cp;
alpha = kTherm ./ rhoCp;                                                   % thermal diffusivity [m^2 / s]


t.start = 0;                                                             % OpenFOAM start time
t.dt = 0.25;                                                              % time step for OpenFOAM controlDict 
t.end = 10;                                              % OpenFOAM end time
t.range = t.start:t.dt:t.end;                                              % range of time

t.plot = t.end;                                                            % specify the time at which the EnKF prediciton is plotted                     
if t.plot > t.end
   error('t.plot (plot time) cannot exceed t.end (end time)!')
end

saveDirectory = '~/dataAssimilation/tests';
EnKF.solverName = 'myLaplacianFoam';                          % OpenFOAM solver name
% EnKF.caseFolder_OF = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF'; % path to the OpenFOAM case directory
EnKF.caseFolder_OF = '~/OpenFOAM/jns14008-5.x/1D_heatConductionEnKF'; 
EnKF.baseLine = '~/OpenFOAM/jns14008-5.x/nominal/postProcessing/sampleDict'; % path to OpenFOAM case directory with steady base line data for plots
t.baseLineData = 1.93;                                                       % time folder of base line data

%     % DecomposePar, only if needed
%     cd(caseFolder_OF);
%     str = sprintf('decomposePar -force -time "0:%.15g"',t.start);   
%     [x,y] = unix(str);

EnKF.varName = cellstr( char('T') );                                         % variables names as they appear in OpenFOAM case time folders
EnKF.N = 100;                                                              % number of cells in the mesh
EnKF.q = 15;                                                               % ensemble size
EnKF.solverRuns = 5;                                                       % Assimilation frequency factor ---> determines number of time steps before measuremetns are assimilated
EnKF.Lradius = 0.1;                                                        % Localization radius [m]

L = 1;
dx = L./EnKF.N;                                                                  % cell spacing
dx2 = dx.*dx;

stdDev.sample = [0];                                                       % standard deviation of initlal sample error of the ensemble
stdDev.w = [1];                                                            % standard deviation for the model
stdDev.v = [1];                                                            % standard deviation for measurements                                        


% Interpolation from states to observations
load('C.mat');
C(5,:)=[];
EnKF.C{1} = C;                                                              % diagonal identity matrices for 1-to-1 mapping for Fv and T

% load('H.mat');
% H(7:9,:)=[];
% EnKF.H{1} = H;
EnKF.H{1} = C;


% Measurements
y.T = [350;400;500;410;500;470;390;330];
% y.T = 500.*ones(8,1); 
EnKF.obsCells{1} = [10 20 30 40 60 70 80 90];
% EnKF.obsCells{1} = 10:10:90; 

EnKF.R{1} = stdDev.v .* eye(numel(y.T));                                    % known measurement covariance matrix

%% === 2.FINITE DIFFERENCE METHOD  === %%

Ti = 300; 
T_FD = Ti.*ones(numel(t.range),EnKF.N+1);   % each row is new time step, each column is new cell edge
T_FD(:,1) = Ti;        % boundary condition
T_FD(:,end) = Ti;      % boundary condition
T_FD(1,:) = Ti;        % Initial condition

B = zeros(EnKF.N+1,1);    % souce % solverRuns = 1;

%input mapping to cells 0.33L and 0.67L
% B( 30:40,1 ) = dt/rhoCp ;
% B( 60:70,1 ) = dt/rhoCp ;
B( round(EnKF.N./3),1 ) = t.dt./rhoCp ;
B( round(EnKF.N.*2./3),1 ) = t.dt./rhoCp ;
uk = zeros(1,numel(t.range));   % tEnd+1 because t=0 stored in column 1, t=tEnd stored in column 1+tEnd

%---------Solve for T Using Runge-Kutta Time Integration---------%
q = 5e7 .* abs(sin(t.range));    % W/m^3
% uk(end+1) = q .* abs(sin(t.end));

for k = 1:numel(t.range)       % TIME LOOP  
    
    for i = 2:EnKF.N     % SPACE LOOP (points are at cell walls, not center)
        
        qSource = B(i) .* q(k);
        
        T_FD(k+1,i) = (t.dt/rhoCp) .* ( (kTherm/dx2) .* (T_FD(k,i-1) - 2.*T_FD(k,i) ...
            + T_FD(k,i+1)) ) + T_FD(k,i) + qSource;
         
    end
    
end
T_FD(end,:) = [];                                                          % eliminate last row which is t.end + t.dt
T_FD = T_FD';                                                              % transpose --> rows = cell index, columns = time

% % PLOT OF TEMPERATURE DISTRIBUTION VS TIME
% fontsize = 17;
% figure();
% imagesc(T_FD');
% xlabel('Cell Number')
% ylabel('Time (s)')
% c = colorbar;
% c.Label.String = 'Temperature (K)';
% % str3 = sprintf('1D Heat Conduction in an Aluminum Rod: q_{source}= 5e+07 * sin(t) W/m^3');
% strTitle = sprintf('1D Heat Conduction in an Aluminum Rod');
% title(strTitle);



% %% Plot contour of T from OpenFOAM simulation myLaplacianFoam
% 
% cd(EnKF.caseFolder_OF);
% 
% nheader = 22;
% 
% T = zeros(EnKF.N,t.end+1);    % rows are each second of time, columns are cell centers
% T(:,1) = 300;       % initialize temperture [K] at t=0, 0 <= x <= L 
% for time = 1:t.end
%     str = sprintf('%g',time);
%     cd(str);
%     fid = fopen('T');
%     temp = textscan(fid,'%f', 'headerlines',nheader);
%     fclose(fid);
%     T(:,time+1) = temp{1};
%     cd ../;
% end
% 
% % PLOT OF TEMPERATURE DISTRIBUTION VS TIME
% fontsize = 17;
% figure();
% imagesc(T');
% xlabel('Cell Number')
% ylabel('Time (s)')
% c = colorbar;
% c.Label.String = 'Temperature (K)';
% % strTitle = sprintf('1D Heat Conduction in an Aluminum Rod: q_{source}= 5e+07 * sin(t) W/m^3');
% strTitle = sprintf('1D Heat Conduction in an Aluminum Rod');
% title(strTitle);


%% === 4. Call EnKF Function ===

% tic                                     
[ filterData,EnKF,t ] = EnKF_Stanford_1DheatConduction_v2( EnKF,stdDev,t,y );  % call ensemble Kalman filter function
% plots( EnKF,stdDev,t,filterData,y,T_FD,saveDirectory )                      % generate plots for temperature along the 1D bar
% toc





%% === TESTING ===s
r = [0.05 0.3, 0.75];

    blue = [0, 0.4470, 0.7410];
    red = [0.6350, 0.0780, 0.1840];
    purple = [0.4940, 0.1840, 0.5560];
    cyan = [0.3010, 0.7450, 0.9330];
    turq = [0, 0.75, 0.75];
    orange = [0.9100,0.4100,0.1700];
    green = [0.4660, 0.6740, 0.1880];
    yellow = [0.9290, 0.6940, 0.1250];

f = figure();
    hold on; box on;
    plot(EnKF.obsCells{1},y.T,'ks','linewidth',1);
    plot(T_FD(:,t.plot ./ t.dt),':','Color',red,'linewidth',1.5);

for i = 1:numel(r)
    EnKF.Lradius = r(i);
    [ filterData,EnKF,t ] = EnKF_Stanford_1DheatConduction_v2( EnKF,stdDev,t,y );  % call ensemble Kalman filter function
    t.kPlot = round( (t.plot - t.start)/(t.dt*EnKF.solverRuns) + 1 );
    Matlab2OF( EnKF,filterData,t.plot,t.kPlot,300 );
   
        plot(filterData{1}(:,t.kPlot),'linewidth',1.5);
    
end
str1 = sprintf('t = %g s, w = %g, v = %g',t.plot,stdDev.w,stdDev.v);
title(str1);
str2 = sprintf('T(t = %g s)', t.plot);
ylabel(str2);
xlabel('Cell ID');
legend('exp.','w/o filter', 'R = 5 cells','R = 30','R = 75'); legend boxoff;
hold off;


saveas(gcf, fullfile(saveDirectory, '1DHeatConductionTEST.tif'));
