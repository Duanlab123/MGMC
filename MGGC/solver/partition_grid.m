function opt=partition_grid(f,alpha,level)
[nr,nc]=size(f);
[M,N]=size(ker(level)); % number of cols in each patch
pr=M;pc=N;
j=level;
npr=floor((nr-1)/2^(j-1)+1);%number of row of patch
npc=floor((nc-1)/2^(j-1)+1);%number of col of patch
row_odd = 1:2:npr;
row_even= 2:2:npr;
col_odd=1:2:npc;
col_even=2:2:npc;
row_odd_ctrx=double(1+(row_odd-1)*2^(j-1));%odd patches location
row_even_ctrx=1+(row_even-1)*2^(j-1);
col_odd_ctrx=1+(col_odd-1)*2^(j-1);
col_even_ctrx=1+(col_even-1)*2^(j-1);
row_odd_idx = zeros(pr*length(row_odd),1);
for ii=1:length(row_odd)
row_odd_idx((ii-1)*pr+1:ii*pr,1)=row_odd_ctrx(ii)-(pr-1)/2:row_odd_ctrx(ii)+(pr-1)/2;
end
row_odd_idx(row_odd_idx<1) = 1;
row_odd_idx(row_odd_idx>nr) = nr;
row_even_idx = zeros(pr*length(row_even),1);
for ii=1:length(row_even)
    row_even_idx((ii-1)*pr+1:ii*pr,1) = row_even_ctrx(ii)-(pr-1)/2:row_even_ctrx(ii)+(pr-1)/2;
end
row_even_idx(row_even_idx<1) = 1;
row_even_idx(row_even_idx>nr) = nr;
col_odd_idx = zeros(pc*length(col_odd),1);
for ii=1:length(col_odd)
    col_odd_idx((ii-1)*pc+1:ii*pc,1) = col_odd_ctrx(ii)-(pc-1)/2:col_odd_ctrx(ii)+(pc-1)/2;
end
col_odd_idx(col_odd_idx<1) = 1;
col_odd_idx(col_odd_idx>nc) = nc;
col_even_idx = zeros(pc*length(col_even),1);
for ii=1:length(col_even)
    col_even_idx((ii-1)*pc+1:ii*pc,1) = col_even_ctrx(ii)-(pc-1)/2:col_even_ctrx(ii)+(pc-1)/2;
end
col_even_idx(col_even_idx<1) = 1;
col_even_idx(col_even_idx>nc) = nc;
%------------------------------------------------------------------------------
pr=pr-2;pc=pc-2;
M = pr; N = pc;
bnd_idx = zeros(4*M, 2);
bnd_idy = zeros(4*M,2);
if M==1 && N==1
    bnd_idy(:,1)=[5;8;6;9];
    bnd_idx(:,1)=bnd_idy(:,1);
    bnd_idx(:,2)=bnd_idx(:,1)-(M+2);
    bnd_idy(:,2)=bnd_idy(:,1)-1;
else
    bnd_idy(1:M+1,1)=(M+2)+2:(M+2):(M+2)+2+(M+2)*(N);
    bnd_idy(M+2:M+2+M,1)=bnd_idy(1:M+1,1)+N;
    bnd_idy(2*M+3:2*M+1+M,1)=(M+2)+3:(M+2)*2-1;
    bnd_idy(3*M+2:4*M,1)=bnd_idy(2*M+3:2*M+1+M,1)+(N)*(M+2);
    bnd_idy(:,2)=bnd_idy(:,1)-1;
    bnd_idx(:,1)=bnd_idy(:,1);
    bnd_idx(:,2)=bnd_idx(:,1)-(M+2);
end
%--------------------------------------------------------------------------------
ntau=pr*pc;
noise=f;
pr=M+2;pc=N+2;
[noise_ch1,noise_ch1_idx,nphi_1,vx_1,vy_1,dv1,s1,mcpara1]=fix_iteration_pro(level,noise,pr,pc,row_odd,col_odd,row_odd_idx,col_odd_idx,bnd_idx,bnd_idy);
[noise_ch2,noise_ch2_idx,nphi_2,vx_2,vy_2,dv2,s2,mcpara2]=fix_iteration_pro(level,noise,pr,pc,row_odd,col_even,row_odd_idx,col_even_idx,bnd_idx,bnd_idy);
[noise_ch3,noise_ch3_idx,nphi_3,vx_3,vy_3,dv3,s3,mcpara3]=fix_iteration_pro(level,noise,pr,pc,row_even,col_odd,row_even_idx,col_odd_idx,bnd_idx,bnd_idy);
[noise_ch4,noise_ch4_idx,nphi_4,vx_4,vy_4,dv4,s4,mcpara4]=fix_iteration_pro(level,noise,pr,pc,row_even,col_even,row_even_idx,col_even_idx,bnd_idx,bnd_idy);
%--------------------------------------------------------------------------------
opt.alpha=alpha;
opt.noise_ch1=noise_ch1;
opt.s1=s1;opt.s2=s2;opt.s3=s3;opt.s4=s4;
opt.nphi_1=nphi_1;
opt.vx_1=vx_1;opt.vy_1=vy_1;opt.dv1=dv1;
opt.noise_ch2=noise_ch2;
opt.nphi_2=nphi_2;opt.vx_2=vx_2;opt.vy_2=vy_2;opt.dv2=dv2;
opt.noise_ch3=noise_ch3;
opt.nphi_3=nphi_3;opt.vx_3=vx_3;opt.vy_3=vy_3;opt.dv3=dv3;
opt.noise_ch4=noise_ch4;
opt.nphi_4=nphi_4;opt.vx_4=vx_4;opt.vy_4=vy_4;opt.dv4=dv4;
opt.ntau=ntau;
opt.npr=npr;opt.npc=npc;
opt.noise=f;
opt.row_odd_idx = row_odd_idx;opt.row_even_idx = row_even_idx;
opt.col_odd_idx = col_odd_idx;opt.col_even_idx = col_even_idx;

opt.row_odd = row_odd;opt.row_even = row_even;
opt.col_odd = col_odd;opt.col_even = col_even;
 
opt.bnd_idx = bnd_idx;
opt.bnd_idy = bnd_idy;

opt.noise_ch1_idx=noise_ch1_idx;
opt.noise_ch2_idx=noise_ch2_idx;
opt.noise_ch3_idx=noise_ch3_idx;
opt.noise_ch4_idx=noise_ch4_idx;
opt.mcpara1=mcpara1;
opt.mcpara2=mcpara2;
opt.mcpara3=mcpara3;
opt.mcpara4=mcpara4;
end


function [noise_ch1,noise_ch1_idx,nphi_1,vx_1,vy_1,dv1,s1,mcpara]=fix_iteration_pro(level,noise,M,N,row_odd,col_odd,row_odd_idx,col_odd_idx,bnd_idx,bnd_idy)
[nr,nc]=size(noise);
z1=noise(row_odd_idx,col_odd_idx);
noise_ch1=zeros(M*N,size(row_odd,2)*size(col_odd,2));
nr=size(row_odd,2);
nc=size(col_odd,2);
for ii = 1:nr
    noise_ch1(:, (ii-1)*nc+1:ii*nc)= reshape(z1((ii-1)*M+1:ii*M,:), M*N, []);
end
phi=ker(level);
phitmp=phi(:);
nphi_1=phitmp*ones(1,size(noise_ch1,2));
s1=sum(nphi_1.*nphi_1,1);
vx_1=nphi_1(bnd_idx(:,1),:)-nphi_1(bnd_idx(:,2),:);
vy_1=nphi_1(bnd_idy(:,1),:)-nphi_1(bnd_idy(:,2),:);
dv1=zeros(M*size(row_odd,2),N*size(col_odd,2));
for ii=1:size(row_odd,2)
    tmp=nphi_1(:,(ii-1)*size(col_odd,2)+1:ii*size(col_odd,2));
    dv1((ii-1)*M+1:ii*N,:)=reshape(tmp,M,[]);      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pr,~]=size(ker(level));
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
x1=ceil(cor1./pr);
y1=cor1-pr.*(x1-1);
x2=ceil(cor2./pr);
y2=cor2-pr.*(x2-1);
x3=ceil(cor3./pr);
y3=cor3-pr.*(x3-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1n=[x1(1:end/2);x2(end/2+1:end)];
y1n=[y1(1:end/2);y2(end/2+1:end)];
x2n=[x2(1:end/2);x1(end/2+1:end)];
y2n=[y2(1:end/2);y1(end/2+1:end)];
x1=x1n;y1=y1n;
x2=x2n;y2=y2n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc=ceil(cor_c./pr);
yc=cor_c-pr.*(xc-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=(x1-x3).*(y2-y3)-(y1-y3).*(x2-x3);
C_matrix=C;%C_matrix*c
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_A=ones(1,size(noise_ch1,2));
C1=C;
C1(C1>0)=1;
C1(C1<0)=-1;
C1=C1*T_A;
C=C*T_A;

%%%%%%%%%%%%%%%%%%%%%%%
para1a=(y1-y3);
para2a=(y2-y3);
para1A=[para1a(1:end/2);-para2a(end/2+1:end)];
para2A=[para2a(1:end/2);-para1a(end/2+1:end)];
%%%%%%%%%%%%%%%%%%%%
para1=para1A*T_A;
para2=para2A*T_A;
%%%%%%%%%%%%%%%%%%%%%%%%%
para1b=(x2-x3);
para2b=(x1-x3);
para1B=[para1b(1:end/2);-para2b(end/2+1:end)];
para2B=[para2b(1:end/2);-para1b(end/2+1:end)];
para3=para1B*T_A;
para4=para2B*T_A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
para1_first=(y1-y3).*(x3-xc)-(x1-x3).*(y3-yc);
para2_first=(y2-y3).*(x3-xc)-(x2-x3).*(y3-yc);
para1_FIST=[para1_first(1:end/2);-para2_first(end/2+1:end)];
para2_FIST=[para2_first(1:end/2);-para1_first(end/2+1:end)];
para5=para1_FIST*T_A;
para6=para2_FIST*T_A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para7=(x3-xc).^2+(y3-yc).^2;
A_part3=para7*T_A;
%%%%%%%%%%%%%%%%%%
mcpara.para1=para1;
mcpara.para2=para2;
mcpara.para3=para3;
mcpara.para4=para4;
mcpara.para5=para5;
mcpara.para6=para6;
mcpara.para7=para7;
mcpara.C=C;
mcpara.C1=0.1*C1;
mcpara.C_matrix=C_matrix;
mcpara.A_part3=A_part3;
[pr,~]=size(ker(level));
idx=pr*((pr-1)/2)+(pr+1)/2;
noise_ch1_idx=noise_ch1(idx,:);

end
