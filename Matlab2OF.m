%% Script to take data from MATLAB and convert to OpenFOAM
% Xinyu Zhao 2017, Joseph Squeo 2019


function Matlab2OF( EnKF,x,time,colIndex,T_BC )

% Write data back into OF format
str = sprintf('%s/%g',EnKF.caseFolder_OF,time);                       % change to time folder of time to be plotted        
cd(str);

for i = 1:EnKF.numVars
    fid = fopen(EnKF.varName{i},'w'); % write permission
    fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\\n');
    fprintf(fid,'| =========                 |                                                 |\n');
    fprintf(fid,'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
    fprintf(fid,'|  \\\\    /   O peration     | Version:  5.x                                   |\n');
    fprintf(fid,'|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n');
    fprintf(fid,'|    \\\\/     M anipulation  |                                                 |\n');
    fprintf(fid,'\\*---------------------------------------------------------------------------*/\n');
    fprintf(fid,'FoamFile\n');
    fprintf(fid,'{\n');
    fprintf(fid,'    version     2.0;\n');
    fprintf(fid,'    format      ascii;\n');
    if strcmp(EnKF.varName{i},'U')
        fprintf(fid,'    class       volVectorField;\n');
    else
        fprintf(fid,'    class       volScalarField;\n');
    end
    fprintf(fid,'    location    "%.15g";\n',time);
    fprintf(fid,'    object       %s;\n',EnKF.varName{i});
    fprintf(fid,'}\n');
    fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
    fprintf(fid,'\n');
    if strcmp(EnKF.varName{i},'T')
        fprintf(fid,'dimensions      [0 0 0 1 0 0 0];\n');
        fprintf(fid,'\n');
        fprintf(fid,'internalField   nonuniform List<scalar>\n');
        fprintf(fid,'%d\n',EnKF.N);
        fprintf(fid,'(\n');
        fprintf(fid,'%.2f\n',x{i}(:,colIndex));
        fprintf(fid,')\n;\n');
        
    elseif strcmp(EnKF.varName{i},'p')
        fprintf(fid,'dimensions      [1 -1 -2 0 0 0 0];\n');
        fprintf(fid,'\n');
        fprintf(fid,'internalField   nonuniform List<scalar>\n');
        fprintf(fid,'%d\n',EnKF.N);
        fprintf(fid,'(\n');
        fprintf(fid,'%.2f\n',x{i}(:,colIndex));
        fprintf(fid,')\n;\n');
    
    else
        fprintf(fid,'dimensions      [0 0 0 0 0 0 0];\n');
        fprintf(fid,'\n');
        fprintf(fid,'internalField   nonuniform List<scalar>\n');
        fprintf(fid,'%d\n',EnKF.N);
        fprintf(fid,'(\n');
        fprintf(fid,'%.9e\n',x{i}(:,colIndex));
        fprintf(fid,')\n;\n');
    end
    
    fprintf(fid,'\nboundaryField\n');
    fprintf(fid,'{\n');
    fprintf(fid,'    inlet\n');
    fprintf(fid,'    {\n');
       if strcmp(EnKF.varName{i},'T')
        fprintf(fid,'        type            fixedValue;\n');
        fprintf(fid,'        value           uniform %.15g;\n',T_BC);
       else
            fprintf(fid,'    type           zeroGradient;\n');
       end
    fprintf(fid,'    }\n');
    fprintf(fid,'    exit\n');
    fprintf(fid,'    {\n');
    if strcmp(EnKF.varName{i},'T')
        fprintf(fid,'        type            fixedValue;\n');
        fprintf(fid,'        value           uniform %.15g;\n',T_BC);
       else
            fprintf(fid,'    type           zeroGradient;\n');
    end
    fprintf(fid,'    }\n');
     fprintf(fid,'    frontAndBack\n');
    fprintf(fid,'    {\n');
    if strcmp(EnKF.varName{i},'T')
        fprintf(fid,'        type            empty;\n');
        %fprintf(fid,'        value           uniform %.15g;\n',T0);
       else
            fprintf(fid,'    type           zeroGradient;\n');
    end
    fprintf(fid,'    }\n');
     fprintf(fid,'    topAndBottom\n');
    fprintf(fid,'    {\n');
    if strcmp(EnKF.varName{i},'T')
        fprintf(fid,'        type            empty;\n');
        %fprintf(fid,'        value           uniform %.15g;\n',T0);
       else
            fprintf(fid,'    type           zeroGradient;\n');
    end
    fprintf(fid,'    }\n');
    fprintf(fid,'}\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'// ************************************************************************* //\n');
    
    fclose(fid);  
end
