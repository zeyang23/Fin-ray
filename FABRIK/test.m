% 21-02-22 ���� �½�FABRIK˼·
% forward and backward reaching inverse kinematics FABRIK
% �����Ǿ���ѧ��ģ�ͺʹ��˶�ѧģ�����ϣ������ڶ�����˿γ���ѧ����˲ʱ���Ե��ص�
% ���Ǿ���ѧ��ģ�ͼ��㸴�Ӷ�̫�ߣ�û��Ҫÿ���㶼�ÿ��Ǿ���ѧ��ģ����������˶�ѧ������
% ��ĳ����̬������һ��С�������ڣ�����˲ʱ���Եĸ����ʹ����������ˣ����˶�ѧҲ��ֻ��Ψһ���

% Extending FABRIK with model constraints ��ƪ���¶����˶�ѧ�Ľ��ܷǳ����������˶�ѧ����ⷽ����Ϊ6��

% ��У�֮ǰ���ſɱȾ�����α�棬�����Ÿ����ؽ���ƽȨ�ģ�ʵ�������ٶȵ�ʱ����Խ�һ����Լ�����Ż�����

clear
clc

n=100;
L=1;

base=[0,0];
target=[0,0.5*L];

d=L/(n-1)*ones(n-1,1);

p0x=transpose(linspace(0,L,n));
p0y=zeros(n,1);

p0=[p0x,p0y];

p=p0;

tol=1e-8;

k=1;

tic;
while(1)
    if norm(p(end,:)-target)<tol
        disp('success')
        disp(k)
        break;
    elseif k>1000
        error('fail to converge')
    else
        p=forward(p,d,target);
        p=backward(p,d,base);
        
        k=k+1;
    end
end
t=toc;
disp(t)

plot_points(p)