function output = addNoise( stdDev,x,numRows,numCols )
% DESCRIPTION:
% addNoise.m adds the appropriate noise to each ensemble of the specified 
% variable, x. the noise is specified in the structure noise.w as an array.
% Additionally, the function filters out negative values by setting these
% values equal to zero.


    
%     for i = 1:EnKF.numVars
        output = x + normrnd(0,stdDev,[numRows numCols]); % add perturbation/noise to the variable x
%     end
                                                
    output(output<0) = 0;                                                   % FILTER out non-physical negative values  

end