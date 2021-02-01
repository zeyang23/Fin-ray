% Class for fixed constrained Rod
classdef fixedRod < handle
    % 基准平面为xy平面
    properties
        % 材料的几何参数
        length;
        width;
        height;
        cross_section;
        
        % 材料的物理参数
        E; % 杨氏模量
        
        % 微元数量
        n_seg;
        
        % 起点的位置和姿态
        start_pos; % 3*1列向量 x y theta theta取角度制
        
        % 终点的位置和姿态
        end_pos;   % 3*1列向量 x y theta theta取角度制
        
        % 根据起点得到的转换后的终点位姿
        conv_end_pos;
        
        % 转换后的 指数坐标所需参量
        conv_q_all;
        conv_w_all;
        conv_xi_all; % 6*n列向量
        conv_xihat_all; % 4*4n列向量
        conv_theta;  % n*1列向量 设置为列向量是为了方便用牛顿法处理
        
        % 实际的 指数坐标所需参量
        q_all;
        w_all;
        xi_all; % 6*n列向量
        xihat_all; % 4*4n列向量
        theta;  % n*1列向量 设置为列向量是为了方便用牛顿法处理
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