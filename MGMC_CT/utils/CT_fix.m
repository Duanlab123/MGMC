function [dw1,C]=CT_fix(w,w1,bnd_idx,bnd_idy,row_odd,col_odd,level,alpha,vx_1,vy_1,dv1,b_part_1_1,b_part_2_1,C,A_phi1)
[pr,pc]=size(ker(level));
U1ch=zeros(pr*pc,size(row_odd,2)*size(col_odd,2));
nr=size(row_odd,2);
nc=size(col_odd,2);
for ii = 1:nr
    U1ch(:, (ii-1)*nc+1:ii*nc) =reshape(w1((ii-1)*pr+1:ii*pr,:), pr*pc, []);
end

[first,b_sec_c_mc]=pre_mean_filter2(U1ch,w,w1,bnd_idx,bnd_idy,row_odd,col_odd,level,alpha,vx_1,vy_1,dv1,b_part_1_1,b_part_2_1,C,A_phi1);
c = new_lanczos_solve(A_phi1,alpha*first' , b_sec_c_mc ,  1e-5,10,level);
 
 
[m,~]=size(C);   
if m~=0
    C1=cell(m+1,1);
    C1(1:m)=C;
    C1(m+1)={c};
    C=C1;
else
    C={c};
end

if level==1
    dw1=(reshape(c',nc,nr))';
else    
dc1=kron(reshape(c',nc,nr)',ones(pr,pc));
dw1=dc1.*dv1;
end
end