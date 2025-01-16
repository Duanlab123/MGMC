function [w,C]=CT_DDM_MC(w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,alpha,vx_1,vy_1,dv1,b_part_1_1,b_part_2_1,C,A_phi1)
add=zeros(size(w));
w1=w(row_odd_idx,col_odd_idx);
[dw1,C]=CT_fix(w,w1,bnd_idx,bnd_idy,row_odd,col_odd,level,alpha,vx_1,vy_1,dv1,b_part_1_1,b_part_2_1,C,A_phi1);
%========================================================
if level==1
    w(row_odd,col_odd)=w(row_odd,col_odd)+dw1;
    add(row_odd,col_odd) = dw1;
else
    if row_odd_idx(1)==row_odd_idx(2)
        dw1(1,:)=dw1(2,:);
    end
    if row_odd_idx(end)==row_odd_idx(end-1)
        dw1(end,:)=dw1(end-1,:);
    end
    if col_odd_idx(1)==col_odd_idx(2)
        dw1(:,1)=dw1(:,2);
    end
    if col_odd_idx(end)==col_odd_idx(end-1)
        dw1(:,end)=dw1(:,end-1);
    end
    w(row_odd_idx,col_odd_idx) = w1+dw1;
    add(row_odd_idx,col_odd_idx) = dw1;
end

end