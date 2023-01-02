function dU1=MMC_fix(w1,bnd_idx,bnd_idy,row_odd,col_even,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,mcpara1,uc1, u1ch_mid1)
[pr,pc]=size(ker(level));
U1ch=zeros(pr*pc,size(row_odd,2)*size(col_even,2));
nr=size(row_odd,2);
nc=size(col_even,2);

for ii = 1:nr
    U1ch(:, (ii-1)*nc+1:ii*nc) =reshape(w1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
end
noisech=noise_ch1;
zbar=noisech-U1ch;
nphi=nphi_1;
zstar=zbar.*(nphi);
zstar=sum(zstar,1)./s1;


mu=1/level;
[dist,add_2,C_matrix,C1]=pre_mean_filter(uc1,pr,mcpara1, u1ch_mid1);
c2=zeros(1,size(U1ch,2));
c3=c2;

for k=1:10
c_old=c2;
add=C_matrix*c2*(-0.75);
add_dist=2*add./add_2;
add_1=add_dist.*C1 ;
dist=0.5*(dist+add_1) ; 
dist_abs=abs(dist);
[~,ind1]=max(dist_abs,[],1);%x
[~,ind2]=min(dist_abs,[],1);%x
ind3=1:size(dist,2);%y(x,y)
index1=ind1+(ind3-1)*size(dist,1);
index2=ind2+(ind3-1)*size(dist,1);
max_dist=dist(index1);
min_dist=dist(index2);
  dm=0.5*max_dist+0.5*min_dist;
%     dm=sum(dist,1)/size(dist,1);
c11=alpha*dm./(s1+alpha );
c22=(s1.*zstar)./(s1+alpha );
c2=   c11  +mu*c22;
add_c=c2-c_old;
c3=c3+c2 ;
if norm(abs(add_c),'inf')<1e-1
    break
end
end
c3=0.5*c3;
ux=U1ch(bnd_idx(:,1),:)-U1ch(bnd_idx(:,2),:);
uy=U1ch(bnd_idy(:,1),:)-U1ch(bnd_idy(:,2),:);
vx=vx_1;vy=vy_1;
z1=zeros(1,size(ux,2));
A1=vx.^2+vy.^2;
a=ones(size(ux,1),1);
C1=ux.*vx+uy.*vy;
tmp=a*c;
z1x=ux+tmp.*vx;
z1y=uy+tmp.*vy;
B1=sqrt(z1x.^2+z1y.^2+1e-4);
T1=sum(A1./B1,1);
C=C1./B1;
T2=sum(C,1);
z1=0.5*(s1*zstar-alpha*T2)./(alpha*T1+s1);
c =c3+z1;
dc1=kron(reshape(c,nc,nr)',ones(pr,pc));
 
dU1=dc1.*dv1;
end
