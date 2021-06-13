% 区分玻璃罐与易拉罐

clear
clc

F_glass=[0,1.8,2.4,2.8,3.2,3.4,4.2,4.5];
F_tin=[0,0.7,1.1,1.3,2.0,2.5,2.9,3.5];

t=0:1:7;

plot(t,F_glass,'-o')
hold on
plot(t,F_tin,'-*')

xlabel('位移指令','FontSize',12)
ylabel('接触力大小/N','FontSize',12)

title('接触力的变化过程','FontSize',18)
legend('玻璃罐','易拉罐','FontSize',12,'Location','NorthWest')