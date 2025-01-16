function [dist] = MC_dist2(U1_ch,nphi,tmpc,level,uc,pr,mcpara, u1ch_mid,vx,vy,bnd_idx,bnd_idy)
%==================================================
para1=mcpara.para1;
para2=mcpara.para2;
para3=mcpara.para3;
para4=mcpara.para4;
para5=mcpara.para5;
para6=mcpara.para6;
C=mcpara.C;
C1=mcpara.C1;
mid_point=mid_point_product(uc,pr, u1ch_mid);
plane=plane_product(pr);
T3=plane(3,:);
U3=mid_point(T3,:);
U1_half=mid_point;
U2_half=flip(mid_point,1);
A_part1=U2_half-U3;
A_part2=U1_half-U3;
A=A_part1.*para1-A_part2.*para2;
B=A_part2.*para3-A_part1.*para4;
UC=ones(size(C,1),1)*uc;
first=A_part1.*para5-A_part2.*para6+C.*(U3-UC);
sec=sqrt(A.^2+B.^2+C.^2);
add_2=sec; 
dist=first./add_2;
dist=(dist.*C1) ;
dist_abs=abs(dist);
[~,ind1]=max(dist_abs,[],1);%x
[~,ind2]=min(dist_abs,[],1);%x
ind3=1:size(dist,2);%y(x,y)
index1=ind1+(ind3-1)*size(dist,1);
index2=ind2+(ind3-1)*size(dist,1);
max_dist=dist(index1);
min_dist=dist(index2);
dm1=(max_dist+ min_dist)/(level*16);
[pr,~]=size(ker(level));
%===================update c to U1ch=============================
U1_ch_new = U1_ch + nphi.*(ones(size(nphi,1),1)*tmpc );
%=================center point BT2 ============================
BT2=U1_ch_new (pr*((pr-1)/2)+(pr+1)/2,:);
%=================distance===============================
u1 = U1_ch_new(1,:);
u2 = U1_ch_new(2,:);
u3 = U1_ch_new(3,:);
u4 = U1_ch_new(4,:);
u5 = U1_ch_new(6,:);
u6 = U1_ch_new(7,:);
u7 = U1_ch_new(8,:);
u8 = U1_ch_new(9,:);
%==============================================
u01 =(u1+u2+u4+BT2)/4;
u02=(u2+BT2)/2;
u03=(u2+BT2+u3+u5)/4;
u04=(u4+BT2)/2;
u05=(BT2+u5)/2;
u06=(u4+BT2+u7+u6)/4;
u07=(BT2+u7)/2;
u08=(BT2+u5+u8+u7)/4;
if level==1
%=================base======================
w1=0.0625;
dist1(1,:)= w1*(BT2-0.5*u02-0.5*u07)./sqrt(( (0.5*u07+0.5*u02-u04)).^2+(0.5*u07-0.5*u02).^2+1);
dist1(2,:)=w1*(BT2-0.5*u02-0.5*u07)./sqrt(( (0.5*u07+0.5*u02-u05)).^2+(0.5*u07-0.5*u02).^2+1);
dist1(3,:)=w1*(BT2-0.5*u04-0.5*u05)./sqrt(( (0.5*u04+0.5*u05-u02)).^2+(0.5*u04-0.5*u05).^2+1);
dist1(4,:)=w1*(BT2-0.5*u04-0.5*u05)./sqrt(( (0.5*u04+0.5*u05-u07)).^2+(0.5*u04-0.5*u05).^2+1);
dist1(5,:)=  (BT2-u2)./sqrt( (BT2-u2).^2+(u2-u3).^2+1);
dist1(6,:)= (BT2-u4)./sqrt( (BT2-u4).^2+(u4-u6).^2+1);
dist1(7,:)=  (BT2-u4)./sqrt( (BT2-u4).^2+(u4-u1).^2+1);
dist1(8,:)= (BT2-u7)./sqrt( (BT2-u7).^2+(u7-u8).^2+1);
%=================base======================
dist2=zeros(8,size(U1_ch_new,2));
dist2(1,:)=  (2*BT2-u2-u4)./sqrt( (BT2-u2).^2+(BT2-u4).^2+1);
dist2(4,:)= (2*BT2-u4-u7)./sqrt( (BT2-u4).^2+(BT2-u7).^2+1);
dist2(7,:)= (2*BT2-u5-u7)./sqrt( (BT2-u5).^2+(BT2-u7).^2+1);
dist2(8,:)= (2*BT2-u2-u5)./sqrt( (BT2-u2).^2+(BT2-u5).^2+1);
dist2(2,:)= (BT2-u5)./sqrt( (BT2-u5).^2+(u5-u3).^2+1);
dist2(3,:)= (BT2-u7)./sqrt( (BT2-u7).^2+(u7-u6).^2+1);
dist2(5,:)= (BT2-u2)./sqrt( (BT2-u2).^2+(u2-u1).^2+1);
dist2(6,:)= (BT2-u5)./sqrt( (BT2-u5).^2+(u5-u8).^2+1); 
dm2=sum(dist1,1)/4+sum(dist2,1)/4; 
dist = dm1+dm2 ;
else
u1_ch_new_x = U1_ch_new(bnd_idx(:,1),:)-U1_ch_new(bnd_idx(:,2),:);
u1_ch_new_y = U1_ch_new(bnd_idy(:,1),:)-U1_ch_new(bnd_idy(:,2),:);
sec_form = sqrt(u1_ch_new_x.^2+u1_ch_new_y.^2+1);
dist_x = (u1_ch_new_x.*vx)./sec_form;
dist_y = (u1_ch_new_y.*vy)./sec_form;
dm2= sum(dist_x + dist_y,1); 
dist =  0.5*dm1+dm2*0.5  ;

end


%===================================================


 
end