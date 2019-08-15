function [Enc,AntiAlias,Direction,scaling_read]=ReadMethods(fn)
fnmethod=[fn,filesep,'method'];
fnacqp=[fn,filesep,'acqp'];
fid = fopen(fnmethod);
R = fread(fid,'char');
fclose(fid);

Methods=(char(R))';

% Find Encoding Matrix size
[id1]=strfind(Methods,'##$PVM_EncMatrix=');
[idend]=id1+100;
EncMatrixString=Methods(id1:idend);
[startIndex,endIndex]=regexp(EncMatrixString,'\d{2,3}');

Enc(1)=str2num(EncMatrixString(startIndex(1):endIndex(1)));
Enc(2)=str2num(EncMatrixString(startIndex(2):endIndex(2)));
Enc(3)=str2num(EncMatrixString(startIndex(3):endIndex(3)));

% Find AntiAlias
[id1]=strfind(Methods,'##$PVM_AntiAlias=');
[idend]=strfind(Methods,'##$PVM_MaxAntiAlias=');
AntiAliasString=Methods(id1:idend);
[newLineindex]=regexp(AntiAliasString,'\n');
[startIndex,endIndex]=regexp(AntiAliasString,'\d');
AntiAlias=AntiAliasString(startIndex(startIndex>newLineindex(1)));
AntiAlias=[str2num(AntiAlias(1)),str2num(AntiAlias(2)),str2num(AntiAlias(3))];

% Find Readout Orientation 
[id1]=strfind(Methods,'##$PVM_SPackArrReadOrient=( 1 )');
SPackArrReadOrientString=Methods(id1:id1+100);
[startIndex,endIndex]=regexp(SPackArrReadOrientString,'\n'); % find either H_F or F_H
Direction=SPackArrReadOrientString(startIndex(1)+1:startIndex(1)+3); % FIND 3 CHARS AFTER THE NEW LINE

%% Open acqp 

fid = fopen(fnacqp);
R2 = fread(fid,'char');
fclose(fid);
acqp=(char(R2))';

[id1]=strfind(acqp,'##$ACQ_scaling_read=');
ScalingString=acqp(id1:id1+30);
[startIndex,endIndex]=regexp(ScalingString,'='); % find either 1 or -1
FirstDigit=ScalingString(startIndex+1);
if FirstDigit=='-'; 
    scaling_read=-1;
else
    scaling_read=1; 
end



end
