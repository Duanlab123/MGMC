function output=MMC_code(f,u0, alpha,max_level)
noise=f;
[nrow,ncol]=size(noise);
u=zeros(nrow,ncol);
num_iter =10;
num_outer=30;
J0=energy_ROF(u,noise,alpha);
Energy(1)=J0;
 
opt0=cell(max_level,1);%store grid in every level
for level=1:max_level
    opt0{level}=partition_grid(noise,alpha,level);
end
%---------------------------------------------------------------
t1=clock;count=2;
for outer_iter=1:num_outer
u_out_old=u;
disp(['outer iteration' num2str(outer_iter)]);        
for level=max_level:-1:1
     opt=opt0{level};
      for ii=1:num_iter               
            u_old=u;
            u=MGMC_solver(u,f,opt,level);               
            error(count-1)=norm(u-u_old,'fro')/norm(u,'fro'); 
            psnr1(count-1)=psnr(uint8(u),uint8(u0));
            count=count+1;                                            
            if error(count-2)<1e-2
                break
            end            
      end            
  
end

Energy(count,1)=energy_ROF(u,noise,alpha);
Energy_out(outer_iter)=energy_ROF(u,noise,alpha);
error_out(outer_iter)=norm(u-u_out_old,'fro')/norm(u,'fro');        
if error_out(outer_iter)<1e-3
     break
end 
end

t2=clock;
t=etime(t2,t1)
Energy=Energy(1:count-1);
error=error(1:count-2);
error_out=error_out(1:outer_iter);
%----------------------------save results----------------------------------
output.u=u;
output.Energy=Energy;
output.Energy_out=Energy_out;
output.error=error;
output.error_out=error_out;
output.iter=count-2;
output.t=t;
end
