classdef finray < handle
    %FIN_RAY 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        pA; % 左底角的固定位置和角度 [xA;yA;alpha]
        pB; % 右底角的固定位置和角度 [xB;yB;beta]
        
        LA; % 左边的长度
        LB; % 右边的长度
        
        % 离散的微元数量
        nA;
        nB;
        
        psi; % 顶角的角度
        
        
        % 物理参数
        wid;
        thi;
        E;
        
        % 内部刚性约束
        constraint_number;
        
        % 3*m矩阵 m为刚性约束的个数
        constraint_ratio;  % 约束长度+左边比例+右边比例 
        constraint_index; % 约束长度+左边对应关节编号+右边
        
        % 左边外力  
        A_force_number % 左边的外力的个数
        
        % 2*p矩阵 p为外力的个数
        A_force_ratio; % 力的大小（压力为正） + 位置比例
        A_force_index; % 力的大小（压力为正） + 关节编号
        
        % 右边外力
        B_force_number % 右边的外力的个数
        
        % 2*q矩阵 q为外力的个数
        B_force_ratio  % 力的大小（压力为正） + 位置比例
        B_force_index  % 力的大小（压力为正） + 关节编号
        
        
        % 计算导数所需的planar_nR对象，对象数组
        RodA;
        RodB;
        
        A_force_array;
        B_force_array;
        
        A_constraint_array;
        B_constraint_array;
    end
    
    methods
        function obj = finray(inputArg1)
            %FIN_RAY 构造此类的实例
            %   此处显示详细说明
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            
        end
    end
end