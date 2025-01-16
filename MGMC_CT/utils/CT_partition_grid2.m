function opt=CT_partition_grid2(A,w,z,patchsize,alpha,level)
[nr,nc]=size(w);
M = patchsize(1); % number of rows in each patch
N = patchsize(2); % number of cols in each patch
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
%%%%%%%%%%%%%%%%%%%%
pr=pr-2;pc=pc-2;
M = pr; N = pc;
bnd_idx = zeros(4*M, 2);
bnd_idy = zeros(4*M,2);
if M==1 && N==1
    bnd_idy(:,1)=[1;2;4;5];
    bnd_idx(:,1)=bnd_idy(:,1);
    bnd_idx(:,2)=bnd_idx(:,1)+(M+2);
    bnd_idy(:,2)=bnd_idy(:,1)+1;
else
    bnd_idy(1:M+1,1)=1:(M+2):1+1+(M+2)*(N);
    bnd_idy(M+2:M+2+M,1)=bnd_idy(1:M+1,1)+N;
    bnd_idy(2*M+3:2*M+1+M,1)=2:M+2-2;
    bnd_idy(3*M+2:4*M,1)=bnd_idy(2*M+3:2*M+1+M,1)+(N)*(M+2);
    bnd_idy(:,2)=bnd_idy(:,1)+1;
    bnd_idx(:,1)=bnd_idy(:,1);
    bnd_idx(:,2)=bnd_idx(:,1)+(M+2);
end
%%%%%%%%%%%%%%%
opt.npr=npr;
opt.npc=npc;

opt.noise=w;
opt.alpha = alpha;

opt.row_odd_idx = row_odd_idx;
opt.row_even_idx = row_even_idx;
opt.col_odd_idx = col_odd_idx;
opt.col_even_idx = col_even_idx;

opt.row_odd = row_odd;
opt.row_even = row_even;
opt.col_odd = col_odd;
opt.col_even = col_even;

opt.patchsize = patchsize;
opt.bnd_idx = bnd_idx;
opt.bnd_idy = bnd_idy;
%%%%%%%%%%%%%%%%%
pr=M+2;pc=N+2;
M=pr;N=pc;


[ vx_1,vy_1,dv1,S1,A_phi1,mcpara1]=fix_iteration_pro2(A,level,w,M,N,row_odd,col_odd,row_odd_idx,col_odd_idx,bnd_idx,bnd_idy);

b_part_1_1=sparse(A_phi1*(z(:)));
save (['b_part_1level' num2str(level)],'b_part_1_1');
add_1={A_phi1,A_phi1'};
save (['add_level' num2str(level)],'add_1');
clear  add_1

[ vx_2,vy_2,dv2,S2,A_phi2,mcpara2]=fix_iteration_pro2(A,level,w,M,N,row_odd,col_even,row_odd_idx,col_even_idx,bnd_idx,bnd_idy);

b_part_1_2=sparse(A_phi2*(z(:)));
save (['b_part_1level' num2str(level)],'b_part_1_2','-append');
add_2={A_phi2,A_phi2'};
save (['add_level' num2str(level)],'add_2','-append');
clear  add_2

[ vx_3,vy_3,dv3,S3,A_phi3,mcpara3]=fix_iteration_pro2(A,level,w,M,N,row_even,col_odd,row_even_idx,col_odd_idx,bnd_idx,bnd_idy);

b_part_1_3=sparse(A_phi3*(z(:)));
save (['b_part_1level' num2str(level)],'b_part_1_3','-append');
add_3={A_phi3,A_phi3'};
save (['add_level' num2str(level)],'add_3','-append');
clear  add_3

[ vx_4,vy_4,dv4,S4,A_phi4,mcpara4]=fix_iteration_pro2(A,level,w,M,N,row_even,col_even,row_even_idx,col_even_idx,bnd_idx,bnd_idy);

b_part_1_4=sparse(A_phi4*(z(:)));
save (['b_part_1level' num2str(level)],'b_part_1_4','-append');
add_4={A_phi4,A_phi4'};
save (['add_level' num2str(level)],'add_4','-append');
clear  add_4

opt.vx_1=vx_1;
opt.vy_1=vy_1;
opt.dv1=dv1;
opt.vx_2=vx_2;
opt.vy_2=vy_2;
opt.dv2=dv2;
opt.vx_3=vx_3;
opt.vy_3=vy_3;
opt.dv3=dv3;
opt.vx_4=vx_4;
opt.vy_4=vy_4;
opt.dv4=dv4;
opt.mcpara1=mcpara1;
opt.mcpara2=mcpara2;
opt.mcpara3=mcpara3;
opt.mcpara4=mcpara4;
end