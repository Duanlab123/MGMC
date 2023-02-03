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
[dist]=pre_mean_filter(uc1,pr,mcpara1, u1ch_mid1,U1ch,bnd_idx,bnd_idy,vx_1,vy_1,alpha,s1,zstar); 
dist_abs=abs(dist);
[~,ind1]=max(dist_abs,[],1);%x
[~,ind2]=min(dist_abs,[],1);%x
ind3=1:size(dist,2);%y(x,y)
index1=ind1+(ind3-1)*size(dist,1);
index2=ind2+(ind3-1)*size(dist,1);
max_dist=dist(index1);
min_dist=dist(index2);
dm=min_dist;
%     dm=sum(dist,1)/size(dist,1);
c11=alpha*dm./(s1+alpha );
c22=(s1.*zstar)./(s1+alpha );
c2=   c11  +mu*c22;
c=c2;
dc1=kron(reshape(c,nc,nr)',ones(pr,pc));
 
dU1=dc1.*dv1;
end

