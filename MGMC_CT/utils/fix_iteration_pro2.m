function [ vx_1,vy_1,dv1,S1,A_phi1,mcpara]=fix_iteration_pro2(A,level,w,M,N,row_odd,col_odd,row_odd_idx,col_odd_idx,bnd_idx,bnd_idy)
[m,n]=size(w);
row_odd_idx(row_odd_idx<1) = 1;
row_odd_idx(row_odd_idx>m) = m;
col_odd_idx(col_odd_idx<1)=1;
col_odd_idx(col_odd_idx>n)=n;
phi=ker(level);
phitmp=phi(:);
nphi_1=phitmp*ones(1,size(row_odd,2)*size(col_odd,2));
vx_1=nphi_1(bnd_idx(:,1),:)-nphi_1(bnd_idx(:,2),:);
vy_1=nphi_1(bnd_idy(:,1),:)-nphi_1(bnd_idy(:,2),:);

dv1=zeros(M*size(row_odd,2),N*size(col_odd,2));
for ii=1:size(row_odd,2)
    tmp=nphi_1(:,(ii-1)*size(col_odd,2)+1:ii*size(col_odd,2));
    dv1((ii-1)*M+1:ii*N,:)=reshape(tmp,M,[]);  
end
[pr,pc]=size(ker(level));
idx_phi1=cell(length(row_odd)*length(col_odd),1);
ii=1;
for q=1:length(row_odd)
    tmp_row_idx=row_odd_idx((q-1)*(pr)+1:q*(pr),1);
    for p=1:length(col_odd)
        tmp_phi=zeros(m,n);
        tmp_col_idx=col_odd_idx((p-1)*(pc)+1:p*(pc),1); 
        tmp_row_idx1=find(tmp_row_idx==1);tmp_row_idx2=find(tmp_row_idx==m);
        tmp_col_idx1=find(tmp_col_idx==1);tmp_col_idx2=find(tmp_col_idx==n);
        tmp_ker=ker(level);
        if ~isempty(tmp_row_idx1)
            tmp_ker(1,:)= tmp_ker(length(tmp_row_idx1),:);
        end
        if ~isempty(tmp_row_idx2)
            tmp_ker(end,:)= tmp_ker(end-length(tmp_row_idx2)+1,:);
        end
        if ~isempty( tmp_col_idx1)
             tmp_ker(:,1)= tmp_ker(:,length( tmp_col_idx1));
        end
        if ~isempty( tmp_col_idx2)
            tmp_ker(:,end)= tmp_ker(:,end-length( tmp_col_idx2)+1);
        end
        tmp_phi(tmp_row_idx,tmp_col_idx)= tmp_ker;
        tmp_idx=find(tmp_phi(:)>0);
        idx_phi1{ii}=tmp_idx;
        ii=ii+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 A_phi1=sparse(size(A,1),length(row_odd)*length(col_odd));
for jj=1:length(row_odd)*length(col_odd)
   tmp_idx= idx_phi1{jj};
   A_phi1(:,jj)=sum(A(:,tmp_idx'),2);
end

A_phi1=A_phi1';
S1=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%
% A_phi1=0;
% S1=0;

%%%%%%%%%%%%%
[nr,nc]=size(w);
row_odd_idx(row_odd_idx<1) = 1;
row_odd_idx(row_odd_idx>nr) = nr;
col_odd_idx(col_odd_idx<1)=1;
col_odd_idx(col_odd_idx>nc)=nc;
z1=w(row_odd_idx,col_odd_idx);
noise_ch1=zeros(M*N,size(row_odd,2)*size(col_odd,2));
nr=size(row_odd,2);
nc=size(col_odd,2);
for ii = 1:nr
    noise_ch1(:, (ii-1)*nc+1:ii*nc)= reshape(z1((ii-1)*M+1:ii*M,:), M*N, []);
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
% x1n=[x1(1:end/2);x2(end/2+1:end)];
% y1n=[y1(1:end/2);y2(end/2+1:end)];
% x2n=[x2(1:end/2);x1(end/2+1:end)];
% y2n=[y2(1:end/2);y1(end/2+1:end)];
% x1=x1n;y1=y1n;
% x2=x2n;y2=y2n;
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
mcpara.C1=C1;
mcpara.C_matrix=C_matrix;
mcpara.A_part3=A_part3;
%%%%%%%%%%%%%%%%%%%%%%
[pr,~]=size(ker(level));
idx=pr*((pr-1)/2)+(pr+1)/2;
noise_ch1_idx=noise_ch1(idx,:);
%%%%output para1 para2 para3 para4 para5 para6 para7 A_part3  C C1 C_matrix

end