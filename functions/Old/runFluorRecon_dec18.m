%%
clear all; close all; clc; 


scanfolders={'L:\basic\divi\Projects\cosart\fluor\fluor\scans\6phantoms_1231\JS_6phantoms_new2.R61\30',...
'L:\basic\divi\Projects\cosart\fluor\fluor\scans\6phantoms_1231\JS_6phantoms_new2.R61\31'};


R=FluorRecon(scanfolders);

R.P.Seq.Sequential=0
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
R.TestMeasOp;
R=R.PseudoInverse;

%%
R=R.ConjugateGradientRecon



