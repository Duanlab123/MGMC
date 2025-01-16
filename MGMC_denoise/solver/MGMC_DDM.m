function [w]=MGMC_DDM(noise,w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,mcpara1,uc1,u1ch_mid1)
w1=w(row_odd_idx,col_odd_idx);

dw1=MGMC_Iteration(w1,bnd_idx,bnd_idy,row_odd,col_odd,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,mcpara1,uc1,u1ch_mid1);
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