% Script demonstrating code usage:
% - Requires path to AIRTOOLS package

clear all; close all; clc;
%addpath('../AIRTools/') % airtools
%addpath('sourceSRS/'); % source 

%%?????????????? ÏÈ×¢ÊÍ ??????????????%%
% --- Phantom Data and Settings --- %
domainDim = 64;
% K = 4;
% noise_lvl = 0.05;
% rngSeed = 1100; 

% xIm = phantomgallery('grains',domainDim, K,rngSeed);
%    xTarget = xIm(:);
    
% lblTarget = ones(domainDim^2,1);
% classVal = 0:1/(K-1):1;
% for i=1:K
%    lblTarget(single(xTarget) == single(classVal(i))) = i; 
% end
%%?????????????? ÏÈ×¢ÊÍ ??????????????%%

% --- Scanner Set-Up --- %
objectSize = 10;
    pixelDim = objectSize/domainDim;
    domainAxis = 0:pixelDim:pixelDim*domainDim;
angProj = 0:1:179;
    Nproj = length(angProj);
dectWidth = objectSize*sqrt(2); %mm
dectNumb  = ceil(domainDim*sqrt(2));
    dectSpace = dectWidth/dectNumb;
    rateUndet = (Nproj*dectNumb)/domainDim^2;
    [A,~,~,~,~,~] = paralleltomo(domainDim,angProj,dectNumb,dectWidth/pixelDim,0);


%%[xRecon,flag,~,returnIter] = cgls(auxA,auxb,0,opt.TolFun,opt.MaxIter,[],mu);
