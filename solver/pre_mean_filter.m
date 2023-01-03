function [dist,add_2,C_matrix,C11]=pre_mean_filter(uc,pr,mcpara,u1ch_mid,U1ch,bnd_idx,bnd_idy,vx_1,vy_1,alpha,s1,zstar)
para1=mcpara.para1;
para2=mcpara.para2;
para3=mcpara.para3;
para4=mcpara.para4;
para5=mcpara.para5;
para6=mcpara.para6;
C11=mcpara.C;
C2=0.1*mcpara.C1;
C_matrix=mcpara.C_matrix;
A_part3=mcpara.A_part3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mid_point=mid_point_product(uc,pr, u1ch_mid);
plane=plane_product(pr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T3=plane(3,:);
U3=mid_point(T3,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U1_half=mid_point;
U2_half=flip(mid_point,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ux=U1ch(bnd_idx(:,1),:)-U1ch(bnd_idx(:,2),:);
uy=U1ch(bnd_idy(:,1),:)-U1ch(bnd_idy(:,2),:);
vx=vx_1;vy=vy_1;
A1=vx.^2+vy.^2;
C1=ux.*vx+uy.*vy;
B1=sqrt(ux.^2+uy.^2+1);
A=A1./B1;
firs=sum(A,1);
C=C1./B1;
sec=sum(C,1);
d1=(-alpha.*sec)./(alpha.*firs+s1);
A_part1=U2_half-U3;
A_part2=U1_half-U3;
A=A_part1.*para1-A_part2.*para2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=A_part2.*para3-A_part1.*para4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UC=ones(size(C11,1),1)*uc;
first=A_part1.*para5-A_part2.*para6+C11.*(U3-UC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sec=sqrt(A.^2+B.^2+C11.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A_part4=(U3-0.5*(U1_half+U2_half)).^2;
% third=A_part3+A_part4;
%%%%%%%%%%%%%%%%%
add_2=sec;%.*third;
dist=first./add_2;
a=ones(size(dist,1),1);
dist=2*(dist.*C2) +a*d1;
%%%%%%%%%%%%%%%%%%%%%%%


end
