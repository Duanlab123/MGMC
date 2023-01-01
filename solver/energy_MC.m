function en = energy_MC(im,alpha,f)
[gx,gy]=mygrad(im);
[gxx,gxy]=mygrad(gx);
[gyx,gyy]=mygrad(gy);
%standard scheme
num = (1+gy.^2).*gxx - gx.*gy.*(gxy+gyx)+ (1+gx.^2).*gyy;
den = (1+gx.^2+gy.^2);
den = sqrt(den).*den*2;
g = num./den;
en1 = sum(abs(g(:)));
data=(f-im).^2;
en=alpha*en1+sum(data(:));
end
function [gx, gy]=mygrad(im)
gx=[im(:,2)-im(:,1) (im(:,3:end)-im(:,1:end-2))./2 im(:,end)-im(:,end-1)];
gy=[im(2,:)-im(1,:) ; (im(3:end,:)-im(1:end-2,:))./2 ; im(end,:)-im(end-1,:)];
end

 