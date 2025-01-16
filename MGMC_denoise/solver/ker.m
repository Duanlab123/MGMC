
%===========================================================
function [k_a ] =ker( level )
if level>1
a=2^(level-1);
tmp=[a/2:-1:0,1:a/2];
 
tmp_x=(1/2).^tmp;
tmp_x=[0,tmp_x,0];
k_a=zeros(size(tmp_x,2),size(tmp_x,2));
for i=1:size(tmp_x,2)
    k_a(i:end-i+1,i:end-i+1)=tmp_x(i);
end
 
else
k_a=[0 0 0; 0 1 0;0 0 0];
end
end

%======================================================
% function [k_a ] =ker( level )
% if level==1
%     k_a=zeros(2^(level-1)+2);
% else
%     k_a=zeros(2^(level-1)+2);
%     [m,n]=size(k_a);
%     k_a=zeros(m+1,n+1);
% end
% k_a(2:end-1,2:end-1)=1;
% end