function sofaStripped = loadSofaFile(hrtf_sofa_path)
% Script using C-style netcdf4 function calls to extract HRIRs and source
% positions from a sofa file

%% Info on sofa file
%ncdisp(hrtf_sofa_path)

%% Open sofa file and determine dimension IDs and lengths
ncid = netcdf.open(hrtf_sofa_path, 'NOWRITE'); % read-only access
for i=1:6% there are 6 dimensions in the sofa standard
    [dimname(i,1), dimlen(i,1)] = netcdf.inqDim(ncid,i-1);
    dimid(i,1) = netcdf.inqDimID(ncid,dimname(i,1));
end  
varname = 'DataType';
DataType = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),varname);

%% Extract IR data
varname = 'Data.IR';
varid = netcdf.inqVarID(ncid,varname); 
[~, ~, dimids, ~] = netcdf.inqVar(ncid,varid);
IR_dim1 = dimlen(dimid(dimids(1,1)+1)+1); %+1s only for MatLab
IR_dim2 = dimlen(dimid(dimids(1,2)+1)+1);
IR_dim3 = dimlen(dimid(dimids(1,3)+1)+1); 
IR = netcdf.getVar(ncid,varid); 
varname = 'Data.SamplingRate';
varid = netcdf.inqVarID(ncid,varname); 
IR_fs = netcdf.getVar(ncid,varid);   

%% Extract positional data
varname = 'SourcePosition';
varid = netcdf.inqVarID(ncid,varname);
[~, ~, dimids, ~] = netcdf.inqVar(ncid,varid);
SourcePosition_dim1 = dimlen(dimid(dimids(1,1)+1)+1); %+1s only for MatLab
SourcePosition_dim2 = dimlen(dimid(dimids(1,2)+1)+1);  
SourcePosition = netcdf.getVar(ncid,varid);

%% Stripped down sofa info
sofaStripped.DataType = DataType;
sofaStripped.IR = IR;
sofaStripped.IR_fs = IR_fs;
sofaStripped.SourcePosition = SourcePosition;
sofaStripped

%% Close sofa file, once info extracted
netcdf.close(ncid);

end