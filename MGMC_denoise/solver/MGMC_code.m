function [ w ,Energy,Energy_out,error,error_out,abs_error] = MGMC_code(z,w,alpha,max_level,u0)
%===================================================
num_outer=2000;
num_iter=10;
Energy= zeros(100*max_level,1);
J_U   = zeros (100*max_level,1);
error = zeros(100*max_level,1);
error_out=zeros(100*max_level,1);
abs_error=zeros(100*max_level,1);
J0=energy_MC(w,z,alpha);
J_U(1)=J0;
Energy(1)=J0;
count=2;iter=2;
%divide the grid
opt0=cell(max_level,1);
for level=1:max_level
    opt0{level}=partition_grid(z,alpha,level);
end
max_loop=1;
for loop=1:max_loop
    disp(['loop number=' num2str(loop)]);
    for outer_iter=1:num_outer
        w_old_out=w;
        disp(['outer iteration' num2str(outer_iter)]);
        for level=max_level:-1:1
            
            disp(['level' num2str(level)]);
            opt=opt0{level};
            J_tmp=zeros(num_iter,1);
            error_tmp=zeros(num_iter,1);
            for ii=1:num_iter
                w_old=w;
                [ w ] =MGMC_MM(w,z,opt,level);
                error(count-1)=norm(w-w_old,'fro')/norm(w,'fro');
                psnr1(count-1)=psnr(uint8(w),uint8(u0));
                error_tmp(ii)=error(count-1);
                count=count+1;
                if error(count-2)<1e-2
                    break
                end
                
            end
            disp(['level=' num2str(level) 'iteration number=' num2str(ii) 'psnr=' num2str( psnr1(count-2))  ]);
            
        end
        J_new =energy_MC(w,z,alpha);
        J_tmp(ii,1)=J_new;
        Energy(count,1)=J_new;
        J_U(iter)=J_new;iter=iter+1;
        error_out(outer_iter)=norm(w-w_old_out,'fro')/norm(w,'fro');
        if error_out(outer_iter)<1e-8
            break
        end
    end
end
Energy=Energy(1:count-1);
error=error(1:count-2);
error_out=error_out(1:outer_iter);
Energy_out=J_U(1:iter-1);
end





