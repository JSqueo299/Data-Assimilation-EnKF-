%% Read OpenFOAM files into Matlab:

function x = OpenFoam2Matlab( EnKF,tNew,j,x )
    
    nheader = 22;                                                          % # of header lines to skip in OF file(s)
    str = sprintf('%s/%.15g',EnKF.caseFolder_OF,tNew);              % change to the new (next time step) time folder
    cd(str);
 
  
    for i = 1:EnKF.numVars
        
        fid = fopen(EnKF.varName{i});                                      % opens the file specified by the variable name
            solverOut(i) = textscan(fid,'%f', 'headerlines',nheader);      % reads OpenFOAM files into MATLAB column vectors
        fclose(fid);
    
        x{i}(:,j) = solverOut{i};
    end
    
end







