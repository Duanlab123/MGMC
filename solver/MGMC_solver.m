function w=MGMC_solver(w,f,opt,level)
alpha = opt.alpha;
row_odd_idx = opt.row_odd_idx;
row_even_idx = opt.row_even_idx;
col_odd_idx = opt.col_odd_idx;
col_even_idx = opt.col_even_idx;
row_odd = opt.row_odd;
row_even = opt.row_even;
col_odd = opt.col_odd;
col_even = opt.col_even;
bnd_idx = opt.bnd_idx;
bnd_idy = opt.bnd_idy;
noise=opt.noise;
s1=opt.s1;
s2=opt.s2;
s3=opt.s3;
s4=opt.s4;
noise_ch1=opt.noise_ch1;nphi_1=opt.nphi_1;vx_1=opt.vx_1;vy_1=opt.vy_1;dv1=opt.dv1;
noise_ch2=opt.noise_ch2;nphi_2=opt.nphi_2;vx_2=opt.vx_2;vy_2=opt.vy_2;dv2=opt.dv2;
noise_ch3=opt.noise_ch3;nphi_3=opt.nphi_3;vx_3=opt.vx_3;vy_3=opt.vy_3;dv3=opt.dv3;
noise_ch4=opt.noise_ch4;nphi_4=opt.nphi_4;vx_4=opt.vx_4;vy_4=opt.vy_4;dv4=opt.dv4;
mcpara1=opt.mcpara1; mcpara2=opt.mcpara2;
mcpara3=opt.mcpara3; mcpara4=opt.mcpara4;
noise_ch1_idx=opt.noise_ch1_idx;noise_ch2_idx=opt.noise_ch2_idx;
noise_ch3_idx=opt.noise_ch3_idx;noise_ch4_idx=opt.noise_ch4_idx;
[uc1,uc2,uc3,uc4,u1ch_mid1,u2ch_mid1,u3ch_mid1,u4ch_mid1]=interp_u1(w,f,level,row_odd,row_even,col_odd,col_even,row_odd_idx,row_even_idx,col_odd_idx,col_even_idx,noise_ch1_idx,noise_ch2_idx,noise_ch3_idx,noise_ch4_idx);
[w]=MMC_DDM_ROF(noise,w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,mcpara1,uc1,u1ch_mid1);
[w]=MMC_DDM_ROF(noise,w,row_odd_idx,col_even_idx,row_odd,col_even,bnd_idx,bnd_idy,level,alpha,noise_ch2,nphi_2,vx_2,vy_2,dv2,s2,mcpara2,uc2,u2ch_mid1);
[w]=MMC_DDM_ROF(noise,w,row_even_idx,col_odd_idx,row_even,col_odd,bnd_idx,bnd_idy,level,alpha,noise_ch3,nphi_3,vx_3,vy_3,dv3,s3,mcpara3,uc3,u3ch_mid1);
[w]=MMC_DDM_ROF(noise,w,row_even_idx,col_even_idx,row_even,col_even,bnd_idx,bnd_idy,level,alpha,noise_ch4,nphi_4,vx_4,vy_4,dv4,s4,mcpara4,uc4,u4ch_mid1);
end

function [w]=MMC_DDM_ROF(noise,w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,mcpara1,uc1,u1ch_mid1)

w1=w(row_odd_idx,col_odd_idx);

dw1=MMC_fix(w1,bnd_idx,bnd_idy,row_odd,col_odd,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,mcpara1,uc1,u1ch_mid1);
if row_odd_idx(1)==row_odd_idx(2)
dw1(1,:)=dw1(2,:);
end
if row_odd_idx(end)==row_odd_idx(end-1)
dw1(end,:)=dw1(end-1,:);
end
if col_odd_idx(1)==col_odd_idx(2)
dw1(:,1)=dw1(:,2);
end
if col_odd_idx(end)==col_odd_idx(end-1)
dw1(:,end)=dw1(:,end-1);
end
w(row_odd_idx,col_odd_idx) = w1+dw1;
end