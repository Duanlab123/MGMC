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


mu=1/level;
[dist,add_2,C_matrix,C1]=pre_mean_filter(uc1,pr,mcpara1, u1ch_mid1,zstar,alpha,s1,ux,uy,vx,vy);
c2=zeros(1,size(ux,2));
c3=c2;

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
c =c3;
dc1=kron(reshape(c,nc,nr)',ones(pr,pc));
 
dU1=dc1.*dv1;
end
