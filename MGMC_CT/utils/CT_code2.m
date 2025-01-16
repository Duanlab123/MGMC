function [w,Energy,Energy_out,error,error_out,t] = CT_code2(A,z,w,alpha,max_level,u0)
Energy= zeros(100*max_level,1);J_U   = zeros (100*max_level,1);
error = zeros(100*max_level,1);error_out=zeros(100*max_level,1);
J0=energy_MC_CT(w,z,alpha,A);J_U(1)=J0;Energy(1)=J0;
count=2;iter=2;
opt0=cell(max_level,1);
for level=max_level:-1:1
    patchsize=size(ker(level));
    opt0{level}=CT_partition_grid2(A,w,z,patchsize,alpha,level);
end
num_outer=2000;
num_iter=10;
%==============================================
for level=max_level:-1:1
    load1='b_part_1level';
    load2=num2str(level);
    output=strcat(load1,load2);
    name_string = ['b_pa_1_level' num2str(level) '=load(output)'];
    eval(name_string);
    load1='add_level';
    load2=num2str(level);
    output=strcat(load1,load2);
    name_string = ['add_level' num2str(level) '=load(output)'];
    eval(name_string);
end
%%%%%%%%%%
[m,n]=size(u0);
%==============================================
t=0;
T1=zeros(size(z,1),1);
for outer_iter=1:num_outer
    w_old_out=w;
    disp(['outer iteration' num2str(outer_iter)]);
for level=max_level:-1:1
    disp(['level' num2str(level)]);
    opt=opt0{level};
    error_tmp=zeros(num_iter,1);
    for ii=1:4
        name_string = ['Tb_part_1_' num2str(ii) '=b_pa_1_level' num2str(level) '.b_part_1_' num2str(ii) ';'];
        eval(name_string)
        name_string = ['add_' num2str(ii) '=add_level' num2str(level) '.add_' num2str(ii) ';'];
        eval(name_string)
    end
    for ii=1:num_iter
        w_old=w;
        t1=clock;
        [w,T1]=MMC_CT2(w,opt,level,Tb_part_1_1,Tb_part_1_2,Tb_part_1_3,Tb_part_1_4,add_1,add_2,add_3,add_4,T1); 
        t2=etime(clock,t1);
        t=t+t2; 
        error(count-1)=norm(w-w_old,'fro')/norm(w,'fro');
        psnr1(count-1)=psnr(uint8(255*w),uint8(255*u0));
        error_tmp(ii)=error(count-1);
        count=count+1;
        if error(count-2)<1e-2
            break
        end
    end
   disp(['level=' num2str(level) 'iteration number=' num2str(ii) 'psnr=' num2str( psnr1(count-2))  ]);
end

J_new =energy_MC_CT(w,z,alpha,A);
J_U(iter)=J_new;iter=iter+1;
error_out(outer_iter)=norm(w-w_old_out,'fro')/norm(w,'fro');
if error_out(outer_iter)<1e-10
    break
end
end
 Energy_out=J_U;
end 

