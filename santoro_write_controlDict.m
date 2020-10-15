%% Script to take data from MATLAB and convert to OpenFOAM
%% Xinyu Zhao 2017, Joseph Squeo 2019

% function tFolder = santoro_write_controlDict( k,EnKF,t )
function santoro_write_controlDict( EnKF,t )

% tFolder = t.start + EnKF.solverRuns .* t.dt .* (k-2);
str1 = sprintf('%s/system/',EnKF.caseFolder_OF);
cd(str1);
maxCourant = 0.5;

%% Write controlDict in OF format

    fid = fopen('controlDict','w'); % write permission
    fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\\n');
    fprintf(fid,'| =========                 |                                                 |\n');
    fprintf(fid,'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
    fprintf(fid,'|  \\\\    /   O peration     | Version:  2.3.0                                 |\n');
    fprintf(fid,'|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n');
    fprintf(fid,'|    \\\\/     M anipulation  |                                                 |\n');
    fprintf(fid,'\\*---------------------------------------------------------------------------*/\n');
    fprintf(fid,'FoamFile\n');
    fprintf(fid,'{\n');
    fprintf(fid,'    version     2.0;\n');
    fprintf(fid,'    format      ascii;\n');
    fprintf(fid,'    class       dictionary;\n');   
    fprintf(fid,'    location    "system";\n');
    fprintf(fid,'    object       controlDict;\n');
    fprintf(fid,'}\n');
    fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
    fprintf(fid,'\n');

    fprintf(fid,'application     %s;\n',EnKF.solverName);
    fprintf(fid,'\n');
    fprintf(fid,'startFrom       startTime;\n');
    fprintf(fid,'\n');
    fprintf(fid,'startTime       %.15g;\n',t.now);
    fprintf(fid,'\n');
    if EnKF.solverRuns == 1
        fprintf(fid,'stopAt          writeNow;\n');
    elseif EnKF.solverRuns >= 2
        fprintf(fid,'stopAt          nextWrite;\n');
    else
        error('solverRuns must be a whole, positive integer!');
    end
    fprintf(fid,'\n');
    fprintf(fid,'endTime         %.15g;\n',t.end);
    fprintf(fid,'\n');
    fprintf(fid,'deltaT          %.15g;\n',t.dt);
    fprintf(fid,'\n');
    fprintf(fid,'writeControl    adjustableRunTime;\n');
    fprintf(fid,'\n');
    if EnKF.solverRuns == 1
        fprintf(fid,'writeInterval   %.15g;\n',t.dt);
    elseif EnKF.solverRuns >= 2
        fprintf(fid,'writeInterval   %.15g;\n',t.dt*EnKF.solverRuns);
    else
        error('EnKF.solverRuns must be a whole, positive integer!');
    end
    fprintf(fid,'\n');
    fprintf(fid,'purgeWrite      0;\n');
    fprintf(fid,'\n');
    fprintf(fid,'writeFormat     ascii;\n');
    fprintf(fid,'\n');
    fprintf(fid,'writePrecision  10;\n');
    fprintf(fid,'\n');
    fprintf(fid,'writeCompression off;\n');
    fprintf(fid,'\n');
    fprintf(fid,'timeFormat      general;\n');
    fprintf(fid,'\n');
    fprintf(fid,'timePrecision   8;\n');
    fprintf(fid,'\n');
    fprintf(fid,'runTimeModifiable true;\n');
    fprintf(fid,'\n');
    fprintf(fid,'adjustTimeStep  no;\n');
    fprintf(fid,'\n');
    fprintf(fid,'maxCo           %.15g;\n',maxCourant);
    fprintf(fid,'\n');
    fprintf(fid,'functions\n');
    fprintf(fid,'{\n');
    fprintf(fid,'}\n');
    fprintf(fid,'\n');
    fprintf(fid,'// ************************************************************************* //\n');

    fclose(fid);  
    cd(EnKF.caseFolder_OF);
end
