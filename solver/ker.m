
function [k_a ] =ker( level )
% if level>1
% a=2^(level-1);
% tmp=[a/2:-1:0,1:a/2];
% tmp_x=(1/2).^tmp;
% x=[0,tmp_x,0];
% y=x;
% k_a=x'*y;
% else
% k_a=[0 0 0; 0 1 0;0 0 0];
% end
% k_a(k_a>0)=1;
% 
if level==1
    k_a=zeros(2^(level-1)+2);
else
    k_a=zeros(2^(level-1)+2);
    [m,n]=size(k_a);
    k_a=zeros(m+1,n+1);
end
k_a(2:end-1,2:end-1)=1;
% aa=2^(level-1);
% tmp=[aa/2:-1:0,1:aa/2];
% len=size(tmp,2);
% ker=zeros(len+2);
% ker_1=zeros((len+1)/2+1);
% ker_2=zeros((len-1)/2+1);
% % ker_3=zeros((len-1)/2+1);
% % ker_4=zeros((len+1)/2+1);
% choose=[1:(len+1)/2]/(len/2+0.5);
% for i=2:(len+1)/2+1
%     tmp_matrix=ones((len+1)/2+1-i+1);
%     tmp_matrix=tmp_matrix*choose(i-1);
%     ker_1(i:end,i:end)=tmp_matrix;
% end
% ker_4=rot90(ker_1,2);
% choose_2=choose(1:end-2);
%     j=1;
% for i=size(choose_2,2):-1:1
%     tmp=choose_2(1:i);
%     tmp_1=zeros(1,(len-1)/2+1-i);
%     tmp_2=[tmp_1,tmp];
%     ker_2(j,:)=tmp_2;
%     j=j+1;   
% end
% ker_3=rot90(ker_2,2);
% ker(1:(len+1)/2+1,1:(len+1)/2+1)=ker_1;
% ker((len+1)/2+1+1:end,1:(len+1)/2)=ker_2;
% ker(1:(len+1)/2,(len+1)/2+2:end)=ker_3;
% ker((len+1)/2+1:end,(len+1)/2+1:end)=ker_4;
% k_a=ker;
% k_a(k_a>0)=1;

end


