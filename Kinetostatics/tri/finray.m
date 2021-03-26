classdef finray < handle
    %FIN_RAY �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
        pA; % ��׽ǵĹ̶�λ�úͽǶ� [xA;yA;alpha]
        pB; % �ҵ׽ǵĹ̶�λ�úͽǶ� [xB;yB;beta]
        
        LA; % ��ߵĳ���
        LB; % �ұߵĳ���
        
        % ��ɢ��΢Ԫ����
        nA;
        nB;
        
        psi; % ���ǵĽǶ�
        
        
        % �������
        wid;
        thi;
        E;
        
        % �ڲ�����Լ��
        constraint_number;
        
        % 3*m���� mΪ����Լ���ĸ���
        constraint_ratio;  % Լ������+��߱���+�ұ߱��� 
        constraint_index; % Լ������+��߶�Ӧ�ؽڱ��+�ұ�
        
        % �������  
        A_force_number % ��ߵ������ĸ���
        
        % 2*p���� pΪ�����ĸ���
        A_force_ratio; % ���Ĵ�С��ѹ��Ϊ���� + λ�ñ���
        A_force_index; % ���Ĵ�С��ѹ��Ϊ���� + �ؽڱ��
        
        % �ұ�����
        B_force_number % �ұߵ������ĸ���
        
        % 2*q���� qΪ�����ĸ���
        B_force_ratio  % ���Ĵ�С��ѹ��Ϊ���� + λ�ñ���
        B_force_index  % ���Ĵ�С��ѹ��Ϊ���� + �ؽڱ��
        
        
        % ���㵼�������planar_nR���󣬶�������
        RodA;
        RodB;
        
        A_force_array;
        B_force_array;
        
        A_constraint_array;
        B_constraint_array;
    end
    
    methods
        function obj = finray(inputArg1)
            %FIN_RAY ��������ʵ��
            %   �˴���ʾ��ϸ˵��
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 �˴���ʾ�йش˷�����ժҪ
            %   �˴���ʾ��ϸ˵��
            
        end
    end
end