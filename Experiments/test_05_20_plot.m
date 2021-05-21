% 处理计算力与测量力的数据并画图

%%
% 有刚性约束 r=20e-3
% delta 5~30

clear
clc

Fsensor_series(:,1)=[0.36;1.45];
Fsensor_series(:,2)=[0.39;2.92];
Fsensor_series(:,3)=[0.43;4.65];
Fsensor_series(:,4)=[0.85;6.50];
Fsensor_series(:,5)=[1.03;9.36];
Fsensor_series(:,6)=[1.10;13.58];

Fcal_norm_series=[];
Fcal_norm_series(:,1)=1.51;
Fcal_norm_series(:,2)=2.95;
Fcal_norm_series(:,3)=4.54;
Fcal_norm_series(:,4)=6.44;
Fcal_norm_series(:,5)=9.10;
Fcal_norm_series(:,6)=13.58;

for i=1:length(Fsensor_series)
    Fsensor_norm_series(:,i)=norm(Fsensor_series(:,i));
end

x=5:5:30;

plot(x,Fsensor_norm_series,'-o')
hold on
plot(x,Fcal_norm_series,'-o')

xlabel('x/mm')
ylabel('F/N')


%%
% 有刚性约束 r=30e-3 31.5
% delta 5~30

clear
clc

Fsensor_series(:,1)=[0.30;1.15];
Fsensor_series(:,2)=[0.37;2.46];
Fsensor_series(:,3)=[0.30;4.19];
Fsensor_series(:,4)=[0.25;6.34];
Fsensor_series(:,5)=[0.35;9.00];
Fsensor_series(:,6)=[0.80;12.60];

Fcal_norm_series=[];
Fcal_norm_series(:,1)=1.24;
Fcal_norm_series(:,2)=2.67;
Fcal_norm_series(:,3)=4.24;
Fcal_norm_series(:,4)=6.06;
Fcal_norm_series(:,5)=8.33;
Fcal_norm_series(:,6)=12.81;

for i=1:length(Fsensor_series)
    Fsensor_norm_series(:,i)=norm(Fsensor_series(:,i));
end

x=5:5:30;

plot(x,Fsensor_norm_series,'-o')
hold on
plot(x,Fcal_norm_series,'-o')

xlabel('x/mm')
ylabel('F/N')

%%
% 无刚性约束 r=30e-3
% delta 5~40

clear
clc

Fsensor_series(:,1)=[0.58;0.94];
Fsensor_series(:,2)=[0.86;2.19];
Fsensor_series(:,3)=[1.24;3.64];
Fsensor_series(:,4)=[1.56;5.14];
Fsensor_series(:,5)=[2.10;6.40];
Fsensor_series(:,6)=[2.70;7.59];
Fsensor_series(:,7)=[3.28;8.57];
Fsensor_series(:,8)=[4.19;9.09];

Fcal_norm_series=[];
Fcal_norm_series(:,1)=1.12;
Fcal_norm_series(:,2)=2.57;
Fcal_norm_series(:,3)=3.90;
Fcal_norm_series(:,4)=5.54;
Fcal_norm_series(:,5)=6.57;
Fcal_norm_series(:,6)=7.57;
Fcal_norm_series(:,7)=8.40;
Fcal_norm_series(:,8)=9.08;

for i=1:length(Fsensor_series)
    Fsensor_norm_series(:,i)=norm(Fsensor_series(:,i));
end

x=5:5:40;

plot(x,Fsensor_norm_series,'-o')
hold on
plot(x,Fcal_norm_series,'-o')

xlabel('x/mm')
ylabel('F/N')

%%
% 有刚性约束 椭圆a=20e-3 b=30e-3 21.5 31.5
% 0:3:18

clear
clc

Fsensor_series(:,1)=[0.70;-0.24];
Fsensor_series(:,2)=[1.40;-0.42];
Fsensor_series(:,3)=[2.22;-0.42];
Fsensor_series(:,4)=[3.26;-0.21];
Fsensor_series(:,5)=[4.52;0.11];
Fsensor_series(:,6)=[6.04;0.52];


Fcal_norm_series=[];
Fcal_norm_series(:,1)=0.70;
Fcal_norm_series(:,2)=1.51;
Fcal_norm_series(:,3)=2.36;
Fcal_norm_series(:,4)=3.27;
Fcal_norm_series(:,5)=4.41;
Fcal_norm_series(:,6)=5.50;

for i=1:length(Fsensor_series)
    Fsensor_norm_series(:,i)=norm(Fsensor_series(:,i));
end

x=3:3:18;

plot(x,Fsensor_norm_series,'-o')
hold on
plot(x,Fcal_norm_series,'-o')

xlabel('x/mm')
ylabel('F/N')