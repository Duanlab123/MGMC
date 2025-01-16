function [w,T1]=MMC_CT2(w,opt,level,b_part_1_1,b_part_1_2,b_part_1_3,b_part_1_4,add_1,add_2,add_3,add_4,T1)
alpha = opt.alpha;
row_odd_idx = opt.row_odd_idx;
row_even_idx = opt.row_even_idx;
col_odd_idx = opt.col_odd_idx;
col_even_idx = opt.col_even_idx;
row_odd = opt.row_odd;
row_even = opt.row_even;
col_odd = opt.col_odd;
col_even = opt.col_even;
bnd_idx = opt.bnd_idx;
bnd_idy = opt.bnd_idy;

vx_1=opt.vx_1;vy_1=opt.vy_1;dv1=opt.dv1;
vx_2=opt.vx_2;vy_2=opt.vy_2;dv2=opt.dv2;
vx_3=opt.vx_3;vy_3=opt.vy_3;dv3=opt.dv3;
vx_4=opt.vx_4;vy_4=opt.vy_4;dv4=opt.dv4;

C=[];

b_part_2_1=add_1{1}*T1;
[w,C]=CT_DDM_MC(w,row_odd_idx,col_odd_idx,row_odd,col_odd,bnd_idx,bnd_idy,level,alpha,vx_1,vy_1,dv1,b_part_1_1,b_part_2_1,C,add_1);

tmp1=add_1{2}*C{1};
 T2=T1+tmp1;

b_part_2_1=add_2{1}*T2;
[w,C]=CT_DDM_MC(w,row_odd_idx,col_even_idx,row_odd,col_even,bnd_idx,bnd_idy,level,alpha,vx_2,vy_2,dv2,b_part_1_2,b_part_2_1,C,add_2);

tmp2=add_2{2}*C{2};
 T3=T2+tmp2;

b_part_2_1=add_3{1}*T3;
[w,C]=CT_DDM_MC(w,row_even_idx,col_odd_idx,row_even,col_odd,bnd_idx,bnd_idy,level,alpha,vx_3,vy_3,dv3,b_part_1_3,b_part_2_1,C,add_3);

tmp3=add_3{2}*C{3};
 T4=T3+tmp3;

b_part_2_1=add_4{1}*T4;
[w,C]=CT_DDM_MC(w,row_even_idx,col_even_idx,row_even,col_even,bnd_idx,bnd_idy,level,alpha,vx_4,vy_4,dv4,b_part_1_4,b_part_2_1,C,add_4);
tmp4=add_4{2}*C{4};
 T1=T4+tmp4;
 

end
