function [EnKF,y,Cx] = preProcessing(EnKF,y,t)

    
    fid = fopen('Cx');
        Cx = textscan(fid,'%f', 'headerlines',22);                            % y-coordinate (radius) cell centrers [m]
        Cx = cell2mat(Cx);
    fclose(fid);

    EnKF.numVars = numel(EnKF.varName);                                        % # of states
    EnKF.numIter = round( (t.end-t.start)/(t.dt*EnKF.solverRuns) );            % calculate the time iterations
    y = struct2cell(y);
    EnKF.numObs = numel(y{1});
end