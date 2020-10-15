function plots( EnKF,noise,t,filterData,y,T_FD,saveDirectory )

% Post-Processing 
     t.kPlot = round( (t.plot - t.start)/(t.dt*EnKF.solverRuns) + 1 );     % plot time index to access correct matricex column vector in filterData
     Matlab2OF( EnKF,filterData,t.plot,t.kPlot,300 );                          % send EnKF state prediction to OpenFOAM before calling sampleDict





% Plot colors
    blue = [0, 0.4470, 0.7410];
    red = [0.6350, 0.0780, 0.1840];
    purple = [0.4940, 0.1840, 0.5560];
    cyan = [0.3010, 0.7450, 0.9330];
    turq = [0, 0.75, 0.75];
    orange = [0.9100,0.4100,0.1700];
    green = [0.4660, 0.6740, 0.1880];
    yellow = [0.9290, 0.6940, 0.1250];

   
% Create Plot(s)
    f = figure();
    hold on; box on;
    plot(EnKF.obsCells{1},y.T,'ks','linewidth',1);
    plot(T_FD(:,t.plot ./ t.dt),':','Color',red,'linewidth',1.5);
    plot(filterData{1}(:,t.kPlot),'Color',turq,'linewidth',1.5);
    legend('measurements','w/o filter', 'filter'); legend boxoff;
%     str1 = sprintf('q = %g, t = %g s, w = %g, v = %g', EnKF.q,tt.lot,noise.w,noise.v);
    str1 = sprintf('1D Heat Conduction with q = %g, w = %g, v = %g', EnKF.q,noise.w,noise.v);
    title(str1);
    str2 = sprintf('T(t = %g s)', t.plot);
    ylabel(str2);
    xlabel('Cell ID');
%     ylim([200 600]);
    hold off;

% saveas(f, '1DHeatConduction.tif');              % save the picture
saveas(gcf, fullfile(saveDirectory, '1DHeatConduction.tif'));
    

end
