clear all
% load forbild_gen_lhw_36_views.mat
addpath utils\
addpath('AIRTools') % airtools
%==================shep=======================
 u_orig = phantom(128); % lambda=1e-1 ;levelmax=4; psnr = 44.0 ~44.44
[M,N] = size(u_orig);           % The discretization points.
theta = 0:10:179;                  % No. of used angles.
p = round(sqrt(2)*N);          % No. of parallel rays.
lambda=1e-1 ;
levelmax=3;
% ===============forbild==============
% u_orig = double(imread('forbild.png'))/255; % lambda=6e-1 ;levelmax=4; psnr = 36.01~38.72
% [M,N] = size(u_orig);           % The discretization points.
% theta = 0:5:179;  % No. of used angles.
% p = round(sqrt(2)*N);           % No. of parallel rays.
% lambda=6e-1 ;
% levelmax=4;
%%===================================================
fprintf(1,'Creating a parallel-beam tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f\n',...
    [N,theta(1),theta(2)-theta(1),theta(end),p]);
% Create the test problem.
[P] = paralleltomo(N,theta,p);
b_ex = P*u_orig(:);
%add noise
randn('seed',2018);
variance = 0.001*max(b_ex(:));
b_ex = b_ex + variance*randn(size(b_ex));
w=zeros(N,N);

t1=clock;
[w,energy,Energy_out,error,error_out,t3]=CT_code2(P,b_ex,w,lambda,levelmax,u_orig);
t2=clock;
t=etime(t2,t1);
figure(100);
imagesc(w); colormap(gray(256));
title(['PSNR=',num2str(psnr(uint8(w*255),uint8(u_orig*255))),'  SSIM=',num2str(ssim(w,u_orig))]);

imwrite(w,'forbild_MCMG_36.png','png');


 