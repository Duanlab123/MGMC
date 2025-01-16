function dU1=MGMC_Iteration(w1,bnd_idx,bnd_idy,row_odd,col_even,level,alpha,noise_ch1,nphi_1,vx_1,vy_1,dv1,s1,mcpara1,uc1, u1ch_mid1)
%====================================================================================
[pr,pc]=size(ker(level));
itmax=5;tauinner=1e-1;
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
%====================================================================================
vx=vx_1;
vy=vy_1;
c=zeros(1,size(U1ch,2));
for iter=1:itmax
    cold=c;
    dist = MC_dist2(U1ch,nphi,c,level,uc1,pr,mcpara1, u1ch_mid1,vx,vy,bnd_idx,bnd_idy);
    f_x= alpha.* dist + s1.*(c-zstar);
    c=c -0.01*f_x/iter;
    if norm(cold-c,'inf')<tauinner
        break
    end
end
dc1=kron(reshape(c,nc,nr)',ones(pr,pc));
dU1=dc1.*dv1;
end

