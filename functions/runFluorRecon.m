%%
clear all; close all; clc; 
%
scanfolders={'L:\basic\divi\Projects\cosart\fluor\fluor\scans\CC_091117_19F_PBS_Jp1\scans\7\BrukerRAW\7',...
'L:\basic\divi\Projects\cosart\fluor\fluor\scans\CC_091117_19F_PBS_Jp1\scans\8\BrukerRAW\8'}% ,...
% 'L:\basic\divi\Projects\cosart\fluor\fluor\scans\CC_091117_19F_PBS_Jp1\scans\9\BrukerRAW\9',...
% 'L:\basic\divi\Projects\cosart\fluor\fluor\scans\CC_091117_19F_PBS_Jp1\scans\10\BrukerRAW\10'};

scanfolders={...
%     'L:\basic\divi\Projects\cosart\fluor\fluor\scans\Archive 2\CC071318_19FPhantom2_Ol1\scans\9\BrukerRAW\9',...
    'L:\basic\divi\Projects\cosart\fluor\fluor\scans\Archive 2\CC071318_19FPhantom2_Ol1\scans\8\BrukerRAW\8',...
    'L:\basic\divi\Projects\cosart\fluor\fluor\scans\Archive 2\CC071318_19FPhantom2_Ol1\scans\7\BrukerRAW\7'}%,...
%     'L:\basic\divi\Projects\cosart\fluor\fluor\scans\Archive 2\CC071318_19FPhantom2_Ol1\scans\6\BrukerRAW\6'}
R=FluorRecon(scanfolders);


R=R.loadBrukerFiles;

R.P.visualization=1;
R.P.zerofilloption=0;
R.Data.improton=[];
R.P.phaseremoval=0

R=R.UpdateP;

R=R.FixDoubleFOV
R=R.ApplyCheckBoard;
R=R.VectorizeData;
R=R.PhaseRemoval;

R=R.LinRecon;
%%
R=R.CalculateSpectra;
% R.TestMeasOp;
R=R.PseudoInverse;

%%
R=R.ConjugateGradientRecon



