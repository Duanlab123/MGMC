function [first,b_sec_c]=pre_mean_filter2(U1_ch,w,w1,bnd_idx,bnd_idy,row_odd,col_odd,level,alpha,vx_1,vy_1,dv1,b_part_1_1,b_part_2_1,C,A_phi1)
[pr,pc]=size(ker(level));
nr=size(row_odd,2);
nc=size(col_odd,2);
vx=vx_1;
vy=vy_1;
%==========================================%
BT2=U1_ch(pr*((pr-1)/2)+(pr+1)/2,:);
uc=BT2;
if level==1
u1 = U1_ch(1,:);
u2 = U1_ch(2,:);
u3 = U1_ch(3,:);
u4 = U1_ch(4,:);
u5 = U1_ch(6,:);
u6 = U1_ch(7,:);
u7 = U1_ch(8,:);
u8 = U1_ch(9,:);
%==============================================
u01 =(u1+u2+u4+BT2)/4;
u02=(u2+BT2)/2;
u03=(u2+BT2+u3+u5)/4;
u04=(u4+BT2)/2;
u05=(BT2+u5)/2;
u06=(u4+BT2+u7+u6)/4;
u07=(BT2+u7)/2;
u08=(BT2+u5+u8+u7)/4;
w1=0.0625;
%==============================================================
dist1(1,:)= w1*(BT2-0.5*u02-0.5*u07)./sqrt(( (0.5*u07+0.5*u02-u04)).^2+(0.5*u07-0.5*u02).^2+1e-7);
dist1(2,:)=w1*(BT2-0.5*u02-0.5*u07)./sqrt(( (0.5*u07+0.5*u02-u05)).^2+(0.5*u07-0.5*u02).^2+1e-7);
dist1(3,:)=w1*(BT2-0.5*u04-0.5*u05)./sqrt(( (0.5*u04+0.5*u05-u02)).^2+(0.5*u04-0.5*u05).^2+1e-7);
dist1(4,:)=w1*(BT2-0.5*u04-0.5*u05)./sqrt(( (0.5*u04+0.5*u05-u07)).^2+(0.5*u04-0.5*u05).^2+1e-7);
dist1(5,:)=  (BT2-u2)./sqrt( (BT2-u2).^2+(u2-u3).^2+1e-7);
dist1(6,:)= (BT2-u4)./sqrt( (BT2-u4).^2+(u4-u6).^2+1e-7);
dist1(7,:)=  (BT2-u4)./sqrt( (BT2-u4).^2+(u4-u1).^2+1e-7);
dist1(8,:)= (BT2-u7)./sqrt( (BT2-u7).^2+(u7-u8).^2+1e-7);
dist2=zeros(8,size(U1_ch,2));
dist2(1,:)=  (2*BT2-u2-u4)./sqrt( (BT2-u2).^2+(BT2-u4).^2+1e-7);
dist2(4,:)= (2*BT2-u4-u7)./sqrt( (BT2-u4).^2+(BT2-u7).^2+1e-7);
dist2(7,:)= (2*BT2-u5-u7)./sqrt( (BT2-u5).^2+(BT2-u7).^2+1e-7);
dist2(8,:)= (2*BT2-u2-u5)./sqrt( (BT2-u2).^2+(BT2-u5).^2+1e-7);
dist2(2,:)= (BT2-u5)./sqrt( (BT2-u5).^2+(u5-u3).^2+1e-7);
dist2(3,:)= (BT2-u7)./sqrt( (BT2-u7).^2+(u7-u6).^2+1e-7);
dist2(5,:)= (BT2-u2)./sqrt( (BT2-u2).^2+(u2-u1).^2+1e-7);
dist2(6,:)= (BT2-u5)./sqrt( (BT2-u5).^2+(u5-u8).^2+1e-7); 
%==============================================================
first1=1./sqrt((BT2-u2).^2+(u2-u3).^2 ++1e-7)+1./sqrt((BT2-u4).^2+(u4-u6).^2   + +1e-7)+1./sqrt((BT2-u5).^2 + (BT2-u7).^2 + +1e-7);
first2=1./sqrt(( (BT2-u2)).^2+(u7-u6).^2+1e-7)+1./sqrt(( (u5-BT2)).^2+(u8-u5).^2++1e-7)+1./sqrt((BT2-u2).^2 + (BT2-u4).^2 + +1e-7);
first3=1./sqrt((BT2-u2).^2 + (BT2-u5).^2 + +1e-7)+1./sqrt((BT2-u4).^2 + (BT2-u7).^2 + +1e-7);
first_1=(first1+first2+first3)/3;
%==============================================================
X1=circshift(w,-1,2);
grad_x_w=w-X1;
Y1=circshift(w,-1);
grad_y_w=w-Y1;
grad_x_w=[zeros(size(w,1),1),grad_x_w];
grad_x_w=[grad_x_w(1,:);grad_x_w];  
grad_y_w=[zeros(1,size(w,1));grad_y_w];
grad_y_w=[grad_y_w(:,1),grad_y_w]; 
ux=zeros(4,size(row_odd,2)*size(col_odd,2));
uy=zeros(4,size(row_odd,2)*size(col_odd,2));
row_odd_l=row_odd+1;col_odd_l=col_odd+1;
TT1=row_odd_l-1;TT2=col_odd_l-1;
T1=grad_x_w(TT1,TT2);T3=grad_x_w(TT1,col_odd_l);
T2=grad_x_w(row_odd_l,TT2);T4=grad_x_w(row_odd_l,col_odd_l);
T1=T1';T2=T2';T3=T3';T4=T4';
ux(1,:)= T1(:)';ux(2,:)= T2(:)';
ux(3,:)= T3(:)';ux(4,:)= T4(:)';
T1=grad_y_w(TT1,TT2);T3=grad_y_w(TT1,col_odd_l);
T2=grad_y_w(row_odd_l,TT2);T4=grad_y_w(row_odd_l,col_odd_l);
T1=T1';T2=T2';T3=T3';T4=T4';
uy(1,:)= T1(:)';uy(2,:)= T2(:)';
uy(3,:)= T3(:)';uy(4,:)= T4(:)';
A1=vx.^2+vy.^2;
%==============================================================
B1=sqrt((ux.^2+uy.^2++1e-7));
E=A1./B1;
first=0.5*first_1+0.5*sum(E,1);
dm=sum(-dist1,1)/8+sum(-dist2,1)/8;
b_part_2=b_part_1_1-b_part_2_1;
b_sec_c = b_part_2+alpha.*dm';
%==============================================================
else
U1_ch_use=zeros(4*pr-4,size(U1_ch,2));
U1_ch_use(1:pr,:)=U1_ch(1:pr,:);
for i=1:pr-2
   tmp=pr+pr*(i-1)+1;
   tmp2=tmp+pr-1;
   tmp_use=pr+2*(i-1)+1;
   tmp_use2=tmp_use+1;
   U1_ch_use([tmp_use;tmp_use2],:)=U1_ch([tmp;tmp2],:);
end
for ii = 1:nr
    U1ch(:, (ii-1)*nc+1:ii*nc) =reshape(w1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
end
ux=U1ch(bnd_idx(:,1),:)-U1ch(bnd_idx(:,2),:);
uy=U1ch(bnd_idy(:,1),:)-U1ch(bnd_idy(:,2),:);   
A1=(vx.^2+vy.^2);
E1=(ux.*vx+uy.*vy);
b_part_2=b_part_1_1-b_part_2_1;
B1=(ux.^2+uy.^2+1e-15);
B1=sqrt(B1);
E=A1./B1;
F=E1./B1;
sec=sum(F,1);
U1_ch_use(tmp_use2+1:end,:)=U1_ch(tmp2+1:end,:);
plane=plane_product(pr);
idx=zeros(pr*2+(pr-2)*2,1);
idx(1:pr,1)=1:pr;
for i=1:pr-2
   idx(pr+2*i-1:pr+2*i,1)=[i*pr+1,(i+1)*pr]'; 
end
idx(end-pr+1:end,1)=pr*(pr-1)+1:pr*pr;
cor1=idx(plane(1,:));
cor2=idx(plane(2,:));
cor3=idx(plane(3,:));
cor_c=pr*((pr-1)/2)+(pr+1)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1=ceil(cor1./pr);
x1=cor1-pr.*(y1-1);
y2=ceil(cor2./pr);
x2=cor2-pr.*(y2-1);
y3=ceil(cor3./pr);
x3=cor3-pr.*(y3-1);
T_A=ones(1,size(U1_ch,2));
x1=x1*T_A;
y1=y1*T_A;
x2=x2*T_A;
y2=y2*T_A;
x3=x3*T_A;
y3=y3*T_A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc=ceil(cor_c./pr);
yc=cor_c-pr.*(xc-1);%norm N=(A,B,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=((x1-x3).*(y2-y3)-(y1-y3).*(x2-x3));
UC=ones(size(C,1),1)*uc;
T1=plane(1,:);
U1=U1_ch_use(T1,:);
T2=plane(2,:);
U2=U1_ch_use(T2,:);
T3=plane(3,:);
U3=U1_ch_use(T3,:);
A=(y1-y3).*U2+(y3-y2).*U1-(y1-y2).*U3;
B=(x3-x1).*U2+(x2-x3).*U1-(x2-x1).*U3;
A=A./abs(C);
B=B./abs(C);
C=C./abs(C);
C1=C;
C1(C1>0)=1;
C1(C1<0)=-1;
dist=0.5*(A.*(xc-x3)+B.*(yc-y3)+C.*(UC-U3))./sqrt(A.^2+B.^2+1);
dist=(dist.*C1) ;
dist=-dist;
dm= sum(dist,1)/(2*size(dist,1));
first=1./sqrt(A.^2+B.^2+1e-15);
first= 0.5*sum(first,1)+0.5*sum(E,1);
b_sec=b_part_2- 2*alpha.*sec'; 
b_sec_c = 0.5*(b_part_2+0.5*alpha.*dm') +0.5*b_sec;
end
 
end