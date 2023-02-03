function en = energy_GC(im,f,alpha)
[gx,gy]=mygrad(im);
[gxx,gxy]=mygrad(gx);
[gyx,gyy]=mygrad(gy);
%standard scheme
num = gxx.*gyy-gxy.*gyx;
den = 1+gx.^2+gy.^2;
den = den.*den;
g = num./den; 
en1 = sum(abs(g(:)));
data=(f-im).^2;
en=alpha*en1+sum(data(:));
end
function [gx, gy]=mygrad(im)
gx=[im(:,2)-im(:,1) (im(:,3:end)-im(:,1:end-2))./2 im(:,end)-im(:,end-1)];
gy=[im(2,:)-im(1,:) ; (im(3:end,:)-im(1:end-2,:))./2 ; im(end,:)-im(end-1,:)];
end