function mid_point=mid_point_product(uc,pr, u1ch_mid)
mid_point=zeros(pr*2+(pr-2)*2,size(uc,2));
%idx_eight=[1,(pr+1)/2,1+pr*(pr-1)/2,(pr+1)/2+pr*(pr-1)/2]';
idx_first_part=[1,pr,pr+2*(pr-2)+1,pr+2*(pr-2)+pr]';
% step1=U1_ch+circshift(U1_ch,-(pr-1)/2);
% step2=(step1+circshift(step1,-pr*(pr-1)/2))./4;
%mid_point(idx_first_part,:)=step2(idx_eight,:);
%%%%%%%%%%%%%%%%%%
mid_point(idx_first_part,:)=u1ch_mid;
%%%%%%%%%%%%%%%%
idx_second=[(pr+1)/2,pr+2*((pr-2)-1)/2+1,pr+2*((pr-2)-1)/2+2,pr+(pr-2)*2+(pr+1)/2]';
second_part_up=mid_point([idx_first_part(1);idx_first_part(3)],:);
second_part_down=mid_point([idx_first_part(2);idx_first_part(4)],:);
mid_point([idx_second(1);idx_second(4)],:)=(second_part_up+second_part_down)./2;

third_part_left=mid_point([idx_first_part(1),idx_first_part(2)]',:);
third_part_right=mid_point([idx_first_part(3),idx_first_part(4)]',:);
mid_point(idx_second(2:3),:)=(third_part_left+third_part_right)./2;

corner1=mid_point(1,:);
corner2=mid_point((pr+1)/2,:);
corner3=mid_point(pr,:);
corner4=mid_point(pr+2*((pr-2)-1)/2+1,:);
corner5=mid_point(pr+2*((pr-2)-1)/2+2,:);
corner6=mid_point(pr+2*(pr-2)+1,:);
corner7=mid_point(pr+2*(pr-2)+(pr+1)/2,:);
corner8=mid_point(pr+2*(pr-2)+pr,:);
if pr==3
else
    tmp=ones((pr+1)/2-2,1);
    tmp_B=[1:(pr+1)/2-2]';
    mid_point(2:(pr+1)/2-1,:)=tmp*corner1+tmp_B*(corner2-corner1)/((pr+1)/2-1);
    mid_point((pr+1)/2+1:pr-1,:)=tmp*corner2+tmp_B*(corner3-corner2)/((pr+1)/2-1);
    mid_point(pr+2*(pr-2)+2:pr+2*(pr-2)+(pr+1)/2-1,:)=tmp*corner6+tmp_B*(corner7-corner6)/((pr+1)/2-1);
    mid_point(pr+2*(pr-2)+(pr+1)/2+1:pr+2*(pr-2)+pr-1,:)=tmp*corner7+tmp_B*(corner8-corner7)/((pr+1)/2-1);
    mid_point(pr+1:2:pr+2*((pr-2)-1)/2+1-1,:)=tmp*corner1+tmp_B*(corner4-corner1)/((pr+1)/2-1);
    mid_point(pr+2:2:pr+2*((pr-2)-1)/2+1,:)=tmp*corner3+tmp_B*(corner5-corner3)/((pr+1)/2-1);
    mid_point(pr+2*((pr-2)-1)/2+3:2:pr+2*(pr-2),:)=tmp*corner4+tmp_B*(corner6-corner4)/((pr+1)/2-1);
    mid_point(pr+2*((pr-2)-1)/2+4:2:pr+2*(pr-2),:)=tmp*corner5+tmp_B*(corner8-corner5)/((pr+1)/2-1);
end

end