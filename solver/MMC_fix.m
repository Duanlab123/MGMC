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
% phi=ker(level);
% phitmp=phi(:);
% nphi=phitmp*ones(1,size(U1ch,2));
zstar=zbar.*(nphi);
%s1=sum(nphi.*nphi,1);
zstar=sum(zstar,1)./s1;

ux=U1ch(bnd_idx(:,1),:)-U1ch(bnd_idx(:,2),:);
uy=U1ch(bnd_idy(:,1),:)-U1ch(bnd_idy(:,2),:);
% vx=nphi(bnd_idx(:,1),:)-nphi(bnd_idx(:,2),:);
% vy=nphi(bnd_idy(:,1),:)-nphi(bnd_idy(:,2),:);
vx=vx_1;
vy=vy_1;
c=zeros(1,size(ux,2));
%%%%%%%%%%%%%%%%%%%%%%%%%iteration%%%%%%%%%%%%%%%%%%%%%%
A1=vx.^2+vy.^2;
a=ones(size(ux,1),1);
C1=ux.*vx+uy.*vy;
for iter=1:5
    tmpc=a*c;
    ucx=ux+tmpc.*vx;
    ucy=uy+tmpc.*vy;
    B1=ucx.^2+ucy.^2+1e-4;
    B1=sqrt(B1);
    A=A1./B1;
    firs=sum(A,1);
    C=C1./B1;
    sec=sum(C,1);
    % c=c-s1.*(c-zstar)+alpha*(firs+sec);
    c=(s1.*zstar-alpha.*sec)./(alpha.*firs+s1);
end
%--------------------------------------------------
mu=1/level;
[dist,add_2,C_matrix,C1]=pre_mean_filter(uc1,pr,mcpara1, u1ch_mid1);
c2=zeros(1,size(ux,2));
c3=c2;
%idx2=C_matrix>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:1
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
c22=(s1.*zstar)./(s1.*firs+alpha );
c2=   c11  +mu*c22;
add_c=c2-c_old;
c3=c3+c2 ;
if norm(abs(add_c),'inf')<1e-1
    break
end
end
c =0.3*c+0.7*c3;
dc1=kron(reshape(c,nc,nr)',ones(pr,pc));
 
dU1=dc1.*dv1;
end