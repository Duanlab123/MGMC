function plane=plane_product(pr)
plane_part_a=zeros(pr-1,2);
plane_part_b=zeros(pr-1,2);
plane_part_a(:,1)=1:pr-1;
plane_part_a(:,2)=1+pr*2+(pr-2)*2-plane_part_a(:,1);
plane_part_b(1:(pr-1)/2,1)=[pr:2:pr+2*(pr-2-1)/2];
plane_part_b(1+(pr-1)/2:end,1)=[pr+2*(pr-2-1)/2+1:-2:pr+1];
plane_part_b(:,2)=1+pr*2+(pr-2)*2-plane_part_b(:,1);
plane_part_a=plane_part_a';
plane_part_b=plane_part_b';
%plane_part_b=[plane_part_b(:,1),plane_part_b(:,end),plane_part_b(:,2:end-1)];
plane_a_group=[plane_part_a;plane_part_b(1,:)];
plane_c_group=[plane_part_a;plane_part_b(2,:)];
plane_c_group=[plane_c_group(2,:);plane_c_group(1,:);plane_c_group(3,:)];

plane_b_group=[plane_part_b;plane_part_a(1,:)];
plane_b_group=[plane_b_group(:,1),plane_b_group(:,end),plane_b_group(:,2:end-1)];

plane_d_group=[plane_part_b;plane_part_a(2,:)];
plane_d_group=[plane_d_group(:,1),plane_d_group(:,end),plane_d_group(:,2:end-1)];
plane_d_group=[plane_d_group(2,:);plane_d_group(1,:);plane_d_group(3,:)];

plane=[plane_a_group,plane_b_group,plane_c_group,plane_d_group];
[~,I]=sort(plane(1,:));
plane=plane(:,I);
end