% Class for fixed constrained Rod
classdef fixedRod < handle
    % ��׼ƽ��Ϊxyƽ��
    properties
        % ���ϵļ��β���
        length;
        width;
        height;
        cross_section;
        
        % ���ϵ��������
        E; % ����ģ��
        
        % ΢Ԫ����
        n_seg;
        
        % ����λ�ú���̬
        start_pos; % 3*1������ x y theta thetaȡ�Ƕ���
        
        % �յ��λ�ú���̬
        end_pos;   % 3*1������ x y theta thetaȡ�Ƕ���
        
        % �������õ���ת������յ�λ��
        conv_end_pos;
        
        % ת����� ָ�������������
        conv_q_all;
        conv_w_all;
        conv_xi_all; % 6*n������
        conv_xihat_all; % 4*4n������
        conv_theta;  % n*1������ ����Ϊ��������Ϊ�˷�����ţ�ٷ�����
        
        % ʵ�ʵ� ָ�������������
        q_all;
        w_all;
        xi_all; % 6*n������
        xihat_all; % 4*4n������
        theta;  % n*1������ ����Ϊ��������Ϊ�˷�����ţ�ٷ�����
    end
    
    methods
        function item = fixedRod(l,w,h,n,pos1,pos2)
        % calculate necessary properties at initialization
            item.length=l;
            item.width=w;
            item.height=h;
            item.n_seg=n;
            item.start_pos=pos1;
            item.end_pos=pos2;
        end
        
        
    end
end