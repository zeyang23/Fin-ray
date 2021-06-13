% ���ֲ�������������

clear
clc

F_glass=[0,1.8,2.4,2.8,3.2,3.4,4.2,4.5];
F_tin=[0,0.7,1.1,1.3,2.0,2.5,2.9,3.5];

t=0:1:7;

plot(t,F_glass,'-o')
hold on
plot(t,F_tin,'-*')

xlabel('λ��ָ��','FontSize',12)
ylabel('�Ӵ�����С/N','FontSize',12)

title('�Ӵ����ı仯����','FontSize',18)
legend('������','������','FontSize',12,'Location','NorthWest')