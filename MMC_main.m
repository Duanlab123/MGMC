%Fast Multi-grid Methods for Minimizing Curvature Energies
% 
clear all
add path solver
load('lena.mat')
u0=double(im);
 
f = u0+0*randn(size(u0));
max_level=4; alpha=15;
tic
output= MMC_code(f,u0,alpha,max_level);  
toc
psnr_u=psnr(uint8(output.u),uint8(u0))
ssim_u=ssim(uint8(output.u),uint8(u0)) 
figure;imshow(uint8(output.u),[])
figure;plot(output.error_out)