close all;
clear; clc;

n=input( 'please input number of points n=');
pic=imread('real.jpg');

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