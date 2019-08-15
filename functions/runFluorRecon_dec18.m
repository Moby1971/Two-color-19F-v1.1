%%
clear all; close all; clc; 
addpath(genpath('L:\basic\divi\Projects\cosart\fluor\fluor\Class'))
addpath(genpath('L:\basic\divi\Projects\cosart\fluor\fluor\functions'))
addpath(genpath('L:\basic\divi\Projects\cosart\Matlab_Collection\spot-master'))
addpath(genpath('L:\basic\divi\Projects\cosart\fluor\fluor\recon3D'))

%
%scanfolders={'L:\basic\divi\Projects\cosart\fluor\fluor\scans\concentrations_1221_newph\JS_F19_6phantomnew.QW1\19',...
%'L:\basic\divi\Projects\cosart\fluor\fluor\scans\concentrations_1221_newph\JS_F19_6phantomnew.QW1\20'};


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



