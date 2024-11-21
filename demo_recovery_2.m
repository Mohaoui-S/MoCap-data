%% Demo of CP decomposition algorithms for MOCAP compeltion problem

clear all;
close all;
clc;

%rng('default');
%% Path to toolboxes
addpath('mocaptoolbox/mocaptoolbox-1.9.8');
addpath('mocaptoolbox/private');
addpath('MoCapToolboxExtension'); 
addpath('aLgorithms_CP');
addpath('tensor_toolbox-master');
%addpath('tensor_toolbox');

%% Initialization and loading data
dataFile = '02_01.c3d';
original = mcread(['data/' dataFile]);

if any(any(isnan(original.data)))
    fprintf('Note that there are already missing data in your file.\n');
    Markes_missing = find(isnan(original.data));
    original.data(Markes_missing)=0;
end
%% Display original data
p = mcinitanimpar;
%  Connexion = [1 2;2 4;4 3;3 1;3 14;4 14;1 13;2 13;14 35;14 23;13 35;13 23;14 38;24 31;24 27;24 29;24 28;31 29;31 27;27 28;28 29;40 41;41 36;36 37; 37 40;36 30;30 25;25 1;... 
%     25 3;34 33;33 39;39 34;39 32;32 26;36 26;22 4;22 2;23 18;18 12;12 5;5 10;10 15;15 8;8 5;35 20;20 9;9 7;7 6;7 17;17 16;16 6]; %FOR 0201
%p.conn =Connexion;
p.conn = mcautomaticbones2(original);    
% fprintf('Automatic connection...\n');
% p.scrsize = [1900/2 800];
% myfighandle = figure(1);
% mc3dplot(original,p,myfighandle);
% title('Original data');

%% Simulate gaps on your data:
Mmissing = 5;  
Fmissing=10; 
Freq=original.freq;
gapsec = Fmissing/Freq; %duration of gaps
%gapsec = .5; %Or you can fix the duration of gaps

incomplete = original;
indices = zeros(size(incomplete.data));
nframes=incomplete.nFrames;
nmarkers = incomplete.nMarkers;
fps = incomplete.freq;
j=randperm(nmarkers);
j=j(1:Mmissing);
gapframe=round(gapsec*fps);
gaps=ones(nmarkers,1);

ind_perm=[];
id1 = 1+floor(rand(1)*(nframes - gapsec*fps));
id2 = min(nframes,id1+gapframe);
for m=j
 ind_perm{m} = zeros(incomplete.nFrames,1);
  for g=1:gaps(m)
      ind_perm{m}(id1:id2)=true;
      indices(id1:id2,(3*m-2):(3*m))=true;  
  end
  ind_perm{m}=logical(ind_perm{m});
end
indices = logical(indices);
incomplete.data(indices)=nan;
Data=original.data;
%% CP 
opts.Rmax =Inf; 
opts.rho =[0, 0, 0];
opts.eps= 1e-3;  
opts.maxiter= 1500;        
opts.tol= 1e-4;        
opts.tol_cnv=1e-6;
opts.gama= 0; %  opts.gama= .001, 0.01; % for spaseCP and smoothCP
opts.model='Cp'; % 'sparse';  'smooth';

%% SPARSITY
opts.rho     = [.01, .01, .01]; % tune this parameter [.0001, 0.0001, .0001];
opts.gama = .001; % tune this parameter 
opts.model='sparse';

%% SMOOTH
opts.rho     = [.001, .001, .001]; % tune this parameter
opts.gama = .01; %0.001,0.0001 % tune this parameter
opts.model='smooth';
%%
R_sequence = CP_algs(incomplete,Data,opts);
Recv=R_sequence.data;
RMSE = sqrt(mean((Recv(:) - Data(:)).^2));
Av_err = mean(abs(Recv(:) - Data(:)));
fprintf ('RMSE %d, average_error %d \r',RMSE, Av_err)

%%