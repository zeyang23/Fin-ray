close all;
clear; clc;

n=10;
pic=imread('delta30.jpg');

I=imshow(pic);

loc_points=zeros(n,2);

%[x,y]=getpts;
for i=1:1:n
    
hold on;  
[x, y]=ginput(1);

hold on;
plot(x,y,'r.')%将点在其中标记出来
 
loc_points(i,1) = x;
loc_points(i,2) = y;

str=['  X:' num2str(x') ', Y:' num2str(y')];
text(x,y,cellstr(str))

end



loc_points_trans=zeros(size(loc_points));
loc_points_trans(:,1)=loc_points(:,1)-loc_points(1,1);
loc_points_trans(:,2)=-(loc_points(:,2)-loc_points(1,2));

lambda=73.24/(loc_points_trans(2,1)-loc_points_trans(1,1));

temp=loc_points_trans*lambda;
loc_points_real=temp(3:end,:);

load('r30_delta30_cal.mat');

rigid_abs_pos-loc_points_real

save('r30_delta30.mat','loc_points_real');