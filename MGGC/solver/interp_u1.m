function [uc1,uc2,uc3,uc4,u1ch_mid,u2ch_mid,u3ch_mid,u4ch_mid]=interp_u1(w,f,level,row_odd,row_even,col_odd,col_even,row_odd_idx,row_even_idx,col_odd_idx,col_even_idx,noisech1_idx,noisech2_idx,noisech3_idx,noisech4_idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n]=size(w);
[m_ker,n_ker]=size(ker(level));
npr=floor((m-1)/2^(level-1)+1);
npc=floor((n-1)/2^(level-1)+1);
row=1:2:npr;
col=1:2:npc;
row_idx=zeros(length(row)*2,1);
col_idx=zeros(length(col)*2,1);
row_ctrx=1+(row -1)*2^(level-1);
col_ctrx=1+(col -1)*2^(level-1);
for ii=1:length(row)
   row_idx((ii-1)*2+1:ii*2,1)=[row_ctrx(ii);row_ctrx(ii)+(m_ker-1)/2]; 
end
row_idx=[1;row_idx;m];
row_idx(row_idx>m)=m;
for ii=1:length(col)
   col_idx((ii-1)*2+1:ii*2,1)=[col_ctrx(ii);col_ctrx(ii)+(n_ker-1)/2]; 
end
col_idx=[1;col_idx;n];
col_idx(col_idx>n)=n;
%%%%%%%%%%%%%%%%%%%%%%%
w_interp=w(row_idx,col_idx);
[m_inter,n_inter]=size(w_interp);
%%%%%%%%%%%%%%%%%%%%%%%%%
[m_inter1,n_inter1]=meshgrid(1:n_inter,1:m_inter);
[x1,y1]=meshgrid(1:0.5:n_inter,1:0.5:m_inter);
w_interp1=interp2(m_inter1,n_inter1,w_interp,x1,y1,'linear');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

row_interp=2:2:size(w_interp1,1)-1;
col_interp=2:2:size(w_interp1,2)-1;
w_opt=w_interp1(row_interp,col_interp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row_odd_opt=1:size(w_opt,1)-1;
row_even_opt=2:size(w_opt,1);
col_odd_opt=1:size(w_opt,2)-1;
col_even_opt=2:size(w_opt,2);
u1ch=w_opt(row_odd_opt,col_odd_opt);
u2ch=w_opt(row_odd_opt,col_even_opt);
u3ch=w_opt(row_even_opt,col_odd_opt);
u4ch=w_opt(row_even_opt,col_even_opt);


 
%%%%%%%%%%%%%%%%%%%% 
tmp_zstar=f-w;
pr=m_ker;
pc=n_ker;
if m_ker~=3
    w1=tmp_zstar(row_odd_idx,col_odd_idx);
    w2=tmp_zstar(row_odd_idx,col_even_idx);
    w3=tmp_zstar(row_even_idx,col_odd_idx);
    w4=tmp_zstar(row_even_idx,col_even_idx);
    
    nr=size(row_odd,2);
    nc=size(col_odd,2);
for ii = 1:nr
    U1_ch(:, (ii-1)*nc+1:ii*nc) =reshape(w1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
    u1ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u1ch((ii-1)*2+1:ii*2,:), 4,[]);
end
%%%%%row_even
    nr=size(row_odd,2);
    nc=size(col_even,2);
for ii = 1:nr
    U2_ch(:, (ii-1)*nc+1:ii*nc) =reshape(w2((ii-1)*pr+1:ii*pr,:), pr*pc, []);
    u2ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u2ch((ii-1)*2+1:ii*2,:), 4,[]);
end   
    nr=size(row_even,2);
    nc=size(col_odd,2);
for ii = 1:nr
    %U3_ch is a change form of w1 with row_even, col_odd, 
    U3_ch(:, (ii-1)*nc+1:ii*nc) =reshape(w3((ii-1)*pr+1:ii*pr,:), pr*pc, []);
    u3ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u3ch((ii-1)*2+1:ii*2,:), 4,[]);    
end
    nr=size(row_even,2);
    nc=size(col_even,2);
for ii = 1:nr
    U4_ch(:, (ii-1)*nc+1:ii*nc) =reshape(w4((ii-1)*pr+1:ii*pr,:), pr*pc, []);
    u4ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u4ch((ii-1)*2+1:ii*2,:), 4,[]);
end    
idx=pr*((pr-1)/2)+(pr+1)/2;
uc1=noisech1_idx-U1_ch(idx,:);
uc2=noisech2_idx-U2_ch(idx,:);
uc3=noisech3_idx-U3_ch(idx,:);
uc4=noisech4_idx-U4_ch(idx,:);

phi=ker(level);
phitmp=phi(:);

nphi1=phitmp*ones(1,size(U1_ch,2));
s1=sum(nphi1.*nphi1,1);

nphi2=phitmp*ones(1,size(U2_ch,2));
s2=sum(nphi2.*nphi2,1);

nphi3=phitmp*ones(1,size(U3_ch,2));
s3=sum(nphi3.*nphi3,1);

nphi4=phitmp*ones(1,size(U4_ch,2));
s4=sum(nphi4.*nphi4,1);

tmp_zstar1=U1_ch.*nphi1;
zstar1=sum(tmp_zstar1,1)./s1; 

tmp_zstar2=U2_ch.*nphi2;
zstar2=sum(tmp_zstar2,1)./s2; 

tmp_zstar3=U3_ch.*nphi3;
zstar3=sum(tmp_zstar3,1)./s3; 

tmp_zstar4=U4_ch.*nphi4;
zstar4=sum(tmp_zstar4,1)./s4; 

else
 
    zstar1=tmp_zstar(row_odd,col_odd);
    zstar2=tmp_zstar(row_odd,col_even);
    zstar3=tmp_zstar(row_even,col_odd);
    zstar4=tmp_zstar(row_even,col_even);
    
    zstar1=zstar1';
    zstar2=zstar2';
    zstar3=zstar3';
    zstar4=zstar4';
    
    zstar1=zstar1(:)';
    zstar2=zstar2(:)';
    zstar3=zstar3(:)';
    zstar4=zstar4(:)';
    
    uc1=noisech1_idx-zstar1;
    uc2=noisech2_idx-zstar2;
    uc3=noisech3_idx-zstar3;
    uc4=noisech4_idx-zstar4;
    
    s1=ones(1,size(zstar1,2));
    s2=ones(1,size( zstar2,2));
    s3=ones(1,size( zstar3,2));
    s4=ones(1,size( zstar4,2));
    
    nr=size(row_odd,2);
    nc=size(col_odd,2);
for ii = 1:nr
    u1ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u1ch((ii-1)*2+1:ii*2,:), 4,[]);
end

    nr=size(row_odd,2);
    nc=size(col_even,2);
for ii = 1:nr
    u2ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u2ch((ii-1)*2+1:ii*2,:), 4,[]);
end

    nr=size(row_even,2);
    nc=size(col_odd,2);
for ii = 1:nr
    u3ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u3ch((ii-1)*2+1:ii*2,:), 4,[]);
end

    nr=size(row_even,2);
    nc=size(col_even,2);
for ii = 1:nr
    u4ch_mid(:, (ii-1)*nc+1:ii*nc)=reshape(u4ch((ii-1)*2+1:ii*2,:), 4,[]);
end    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end