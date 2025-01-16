clear all
addpath solver\
addpath data\
%load shapes.mat %max_level=2; alpha=20; psnr = 36.05    
 %load synthetic.mat %max_level=2; alpha=35; psnr = 38.44 
%load square.mat %max_level=3; alpha=35; psnr = 40.70 
load peppers.mat %max_level=3; alpha=15; psnr = 32.11 
%load man.mat %max_level=3; alpha=15; psnr = 31.65 
%load flower.mat %max_level=3; alpha=15; psnr = 34.55 
u0 = double(clean*255);
f = double(noisy*255);
lambda =15;
levelmax=3;
w=zeros(size(u0));
t1=clock;
[w,energy,Energy_out,error,error_out]=MGMC_code(f,w,lambda,levelmax,u0);
t2=clock;
t=etime(t2,t1);
psnr_u=psnr(uint8(w),uint8(u0))
 