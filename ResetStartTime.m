
%% Calling OpenFOAM from Matlab:
function ResetStartTime(EnKF,time)

    cd(EnKF.caseFolder_OF);      % change directory to the case foler
    str1 = sprintf('cp T_0Folder T');
    str2 = sprintf('mv T ./%.15g',time);
    [dummy1,dummy2] = unix(str1);
    [dummy1,dummy2] = unix(str2);  
    
end







