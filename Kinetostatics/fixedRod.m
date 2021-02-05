% Class for fixed constrained Rod
% 21-02-04�������⣬��ǰĩ��ֻ������ʼ�����Ϸ�������һ���ޣ�������������ޣ����޷���Ԥ�ϵ�һ����⡣
% ��Ϊԭ�㴦û�йؽڣ���һ���ؽ���x���������ϡ�˵��ͼ��fixedRod_02_04.jpg

% 21-02-05 �������ƺ��������⡣�ſɱȾ����Ƿ��һ������任

classdef fixedRod < handle
    % ��׼ƽ��Ϊxyƽ��
    properties
        % ���ϵļ��β���
        Ltotal;
        wid;
        thi;
        cross_section;
        Iz; % ת������
        
        % ���ϵ��������
        E; % ����ģ��
        
        % ΢Ԫ�����볤��
        n_seg;
        seg_length;
        
        % �նȾ���
        K_theta;
        
        % ����λ�ú���̬
        start_pos; % 3*1������ x y theta thetaȡ�Ƕ���
        
        % �յ��λ�ú���̬
        end_pos;   % 3*1������ x y theta thetaȡ�Ƕ���
        
        % �������õ���ת������յ�λ��
        conv_end_pos;
        conv_gt;
        
        % ת����� ָ�������������
        conv_q_all;
        conv_w_all;
        conv_xi_all; % 6*n������
        conv_xihat_all; % 4*4n������
        conv_g0;
        
        conv_theta;  % n*1������ ����Ϊ��������Ϊ�˷�����ţ�ٷ����� ��ע���ǻ�����
        conv_F; % ĩ�˵�Լ���� 6*1������
        
        conv_g;
        conv_jacobs;
        conv_K_J;
        
        conv_pos_all % ����n+2��������� 2*(n+2)
        
        % ʵ�ʵ�
        pos_all; % ����n+2��������� 2*(n+2)
    end
    
    methods
       %% ���캯��
        function item = fixedRod(E,l,w,h,n,pos1,pos2)
        % calculate necessary properties at initialization
            item.E=E;
            item.Ltotal=l;
            item.wid=w;
            item.thi=h;
            item.cross_section=item.wid*item.thi;
            item.Iz=1/12*item.thi.^3*item.wid;
            
            item.n_seg=n;
            item.start_pos=pos1;
            item.end_pos=pos2;
            
            item.seg_length=item.Ltotal/item.n_seg;
            
            item.K_theta=diag((item.E*item.Iz/item.seg_length)*ones(n,1));
        end
        
       %% ��ʼ��conv����
        function init_exp(obj)
            theta_start=obj.start_pos(3);
            theta_start_rad=theta_start/180*pi;
            delta_x=obj.end_pos(1)-obj.start_pos(1);
            delta_y=obj.end_pos(2)-obj.start_pos(2);
            obj.conv_end_pos(1)=delta_x*cos(theta_start_rad)+delta_y*sin(theta_start_rad);
            obj.conv_end_pos(2)=delta_y*cos(theta_start_rad)-delta_x*sin(theta_start_rad);
            obj.conv_end_pos(3)=obj.end_pos(3)-obj.start_pos(3);
            
            obj.conv_gt=[rotz(obj.conv_end_pos(3)),[obj.conv_end_pos(1);obj.conv_end_pos(2);0];0,0,0,1];
            
            obj.conv_q_all=zeros(3,obj.n_seg);
            for i=1:obj.n_seg
                obj.conv_q_all(:,i)=[1;0;0]*obj.seg_length*(i-1/2);
            end
            
            obj.conv_w_all=[];
            for i=1:obj.n_seg
                obj.conv_w_all=[obj.conv_w_all,[0;0;1]];
            end
            
            obj.conv_g0=[eye(3),[obj.Ltotal;0;0];0,0,0,1];
            
            obj.conv_xi_all=zeros(6,obj.n_seg);
            obj.conv_xihat_all=zeros(4,4*obj.n_seg);
            for a=1:obj.n_seg
                obj.conv_xi_all(1:3,a)=-obj.cha(obj.conv_w_all(1:3,a),obj.conv_q_all(1:3,a));
                obj.conv_xi_all(4:6,a)=obj.conv_w_all(:,a);
                obj.conv_xihat_all(1:3,4*(a-1)+1:4*(a-1)+3)=[0,-obj.conv_w_all(3,a),obj.conv_w_all(2,a);obj.conv_w_all(3,a),0,-obj.conv_w_all(1,a);-obj.conv_w_all(2,a),obj.conv_w_all(1,a),0];
                obj.conv_xihat_all(1:3,4*a)=-obj.cha(obj.conv_w_all(1:3,a),obj.conv_q_all(1:3,a));
            end
           
            obj.conv_theta=zeros(obj.n_seg,1);
            [obj.conv_jacobs,~]=obj.exp_jacob(obj.conv_w_all,obj.conv_q_all,obj.conv_g0,obj.conv_theta);
            obj.conv_g=obj.exp_fkine(obj.conv_w_all,obj.conv_q_all,obj.conv_g0,obj.conv_theta); 
            
            obj.conv_F=zeros(6,1);
            
            obj.conv_K_J=-obj.partial_J_theta(obj.conv_w_all,obj.conv_q_all,obj.conv_jacobs,obj.conv_F);
        end
        
        
       %% ����conv�µĲ���
       function update_conv(obj)
            [obj.conv_jacobs,~]=obj.exp_jacob(obj.conv_w_all,obj.conv_q_all,obj.conv_g0,obj.conv_theta);
            obj.conv_g=obj.exp_fkine(obj.conv_w_all,obj.conv_q_all,obj.conv_g0,obj.conv_theta);
            obj.conv_K_J=-obj.partial_J_theta(obj.conv_w_all,obj.conv_q_all,obj.conv_jacobs,obj.conv_F);
       end
  
        
       %% ţ�ٷ����
       function Newton_conv(obj,TOL)
            k=1;
            while(1)
                if k>500
                    error("can not converge")
                end
                J=cal_Jacob_conv(obj);
                b=cal_constraint_conv(obj);
                delta=-pinv(J)*b;
                
                % ����conv����ϵ�µĲ���
                obj.conv_theta=obj.conv_theta+delta(1:obj.n_seg);
                obj.conv_F=obj.conv_F+delta(end-5:end);
                
                update_conv(obj);
                
                k=k+1;
        %         if(norm(x(:,k)-x(:,k-1))<TOL)
                if(norm(cal_constraint_conv(obj))<TOL)
                    fprintf('Newton Method converge: iteration = %d\n',k-1)
                    fprintf('norm(e) = %E\n',norm(cal_constraint_conv(obj)))
                    break;
                end
            end
       end
        
       %% ����Լ�� ��conv����ϵ��
       function r=cal_constraint_conv(obj)

%             [jacobs,~]=obj.exp_jacob(obj.conv_w_all,obj.conv_q_all,obj.conv_g0,obj.conv_theta);
% 
%             g=obj.exp_fkine(obj.conv_w_all,obj.conv_q_all,obj.conv_g0,obj.conv_theta); 

            tau=obj.K_theta*obj.conv_theta-transpose(obj.conv_jacobs)*obj.conv_F;

            error=logm(obj.conv_g/obj.conv_gt);
            e(1:3,1)=error(1:3,4);
            e(4,1)=-error(2,3);
            e(5,1)=error(1,3);
            e(6,1)=-error(1,2);

            r=[e;tau];

       end
        
       %% �����ſɱȾ��� ��conv����ϵ��
        function JACOBIANs=cal_Jacob_conv(obj)

            JACOBIANs=zeros(6+obj.n_seg);
%             JACOBIANs(1:6,1:obj.n_seg)=obj.adjoint(inv(obj.conv_gt))*obj.conv_jacobs;
            JACOBIANs(1:6,1:obj.n_seg)=obj.conv_jacobs;
            JACOBIANs(7:end,1:obj.n_seg)=obj.K_theta+obj.conv_K_J;
            JACOBIANs(7:end,obj.n_seg+1:end)=-transpose(obj.conv_jacobs);
        end
        
        function A=partial_J_theta(obj,w_all,q_all,Jacobs,F)
            xi_a=[-obj.cha(w_all,q_all);w_all];
            n=size(xi_a,2);
            A=zeros(n,n);
            for J=2:n
                for K=1:J-1
                    A(J,K)=transpose(obj.Z_adjoint(Jacobs(:,K))*Jacobs(:,J))*F;
                end
            end
        end
        
        
       %% һЩ���ܺ��� exp_jacob exp_fkine
        function [jacob0,jacobe] = exp_jacob(obj,w_all,q_all,g0,theta_all)
            n=length(theta_all);
            jacob0=zeros(6,size(w_all,2));
            xi_hat_all=obj.hight(w_all,q_all);                          %��R^6��se(3)
            xi_a=[-obj.cha(w_all,q_all);w_all];

            mat_A=zeros(4,4*n);
            mat_A(:,1:4)=eye(4);
            for i = 2:n
                mat_A(:,4*i-3:4*i)= mat_A(:,4*i-7:4*i-4)*obj.mult(xi_hat_all(1:4,4*i-7:4*i-4),theta_all(i-1));
            end

            A_end=mat_A(:,4*n-3:4*n)*obj.mult(xi_hat_all(:,4*n-3:4*n),theta_all(n))*g0;

            for a=1:size(w_all,2)
                jacob0(:,a)=obj.adjoint(mat_A(:,4*(a-1)+1:4*a))*xi_a(:,a);
            end

            jacobe=obj.adjoint(A_end)\jacob0;
        end

        function [xi_all] = hight(obj,w_all,q_all)       %%��R^6��se(3)�������ϴ�ֱ��˲���
            xi_all=zeros(4,4*size(w_all,2));
            for a=1:size(w_all,2)
                xi_all(1:3,4*(a-1)+1:4*(a-1)+3)=[0,-w_all(3,a),w_all(2,a);w_all(3,a),0,-w_all(1,a);-w_all(2,a),w_all(1,a),0];
                xi_all(1:3,4*a)=-obj.cha(w_all(1:3,a),q_all(1:3,a));
            end
        end
        
        function [r] = cha(obj,a,b)
            r=zeros(3,size(a,2));
            for i=1:size(a,2)
                hight_a=[0 -a(3,i) a(2,i);a(3,i) 0 -a(1,i);-a(2,i) a(1,i) 0];
                r(:,i)=hight_a*b(:,i);
            end
        end
        
        
        function [gi_mult] = mult(obj,xihat_all,theta_all)
            gi_mult=eye(4);
            for a=1:length(theta_all)
                w=xihat_all(1:3,4*(a-1)+1:4*(a-1)+3);
                v=xihat_all(1:3,4*a);
                exp_w=(eye(3)+w*sin(theta_all(a))+w*w*(1-cos(theta_all(a))));
                x=(eye(3)-exp_w)*w*v+[w(3,2),w(1,3),w(2,1)]*v*transpose([w(3,2),w(1,3),w(2,1)])*theta_all(a);
                gi_mult=gi_mult*[exp_w,x;0 0 0 1];
            end
        end

        function g_st = exp_fkine(obj,w_all,q_all,g0,theta_all)%�����˶�ѧ����
            xihat=obj.hight(w_all,q_all);                          %��R^6��se(3)
            gi_mult=obj.mult(xihat,theta_all);                     %��se(3)��SE(3)���۳�
            g_st=gi_mult*g0;                                    %��ֵ�������
        end
        
        function [adj] = adjoint(obj,A)
            a=A(1:3,4);ph=[0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
            adj=[A(1:3,1:3),ph*A(1:3,1:3);zeros(3),A(1:3,1:3)];
        end
        
        function B=Z_adjoint(obj,xi)
        % 01/26���£�ԭ�ȵ�Z_adjoint�Ǵ�ģ�������[omega;v]��
        % Z_adjoint2����[v;omega]
        % 01/26����,Z_adjoint2Ҳ�Ǵ�ģ�û�к�F��Ӧ��ԭ���е�FΪ[m;f]�����õ���[f;m]
            omega=xi(4:6);
            v=xi(1:3);
            B=zeros(6,6);
            B(1:3,1:3)=obj.xev(omega);
            B(4:6,4:6)=obj.xev(omega);
            B(1:3,4:6)=obj.xev(v);
        end
        
        function vec=vex(obj,mat)
            % inversed skew-symmetric transform
            vec = [mat(6), -mat(3), mat(2)];
        end
        
        function mat=xev(obj,vec)
            % skew-symmetric transform

          mat = [0, -vec(3), vec(2)
                vec(3), 0, -vec(1)
                -vec(2), vec(1), 0];
        end 
        
        
       %% ���ӻ�����
       function plot_pos_conv(obj)
            theta_solve=obj.conv_theta;
            Theta=zeros(1,length(theta_solve));
            for k=1:length(Theta)
                Theta(k)=theta_solve(k);
            end

            g_exp=obj.conv_g;

            %�õ������˵�λ��
            pos=[0;0];
            pos=[pos,[obj.conv_q_all(1,1);obj.conv_q_all(2,1)]];
            for i =2:length(Theta)
                g0i=[eye(3),[(i-0.5)*obj.seg_length;0;0];0,0,0,1];
                gsti=obj.exp_fkine(obj.conv_w_all(:,1:i-1),obj.conv_q_all(:,1:i-1),g0i,Theta(:,1:i-1));
                pos=[pos,gsti(1:2,4)];
            end
            pos=[pos,g_exp(1:2,4)];

            plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)


            new_xy=g_exp(1:2,1:2)*[1 0;0 1];
            new_x=new_xy(1,:);
            new_y=new_xy(2,:);
            hold on
            End_pos=g_exp(1:2,4);
            location_x=[End_pos(1),End_pos(1)];
            location_y=[End_pos(2),End_pos(2)];
            quiver(location_x,location_y,new_x,new_y,0.05*obj.Ltotal)
            axis equal

       end
       
       function plot_pos(obj)
            theta_solve=obj.conv_theta;
            Theta=zeros(1,length(theta_solve));
            for k=1:length(Theta)
                Theta(k)=theta_solve(k);
            end

            g_exp=obj.conv_g;

            %�õ������˵�λ��
            pos=[0;0];
            pos=[pos,[obj.conv_q_all(1,1);obj.conv_q_all(2,1)]];
            for i =2:length(Theta)
                g0i=[eye(3),[(i-0.5)*obj.seg_length;0;0];0,0,0,1];
                gsti=obj.exp_fkine(obj.conv_w_all(:,1:i-1),obj.conv_q_all(:,1:i-1),g0i,Theta(:,1:i-1));
                pos=[pos,gsti(1:2,4)];
            end
            pos=[pos,g_exp(1:2,4)];
            
            for j=1:length(pos)
                R=rotz(obj.start_pos(3));
                pos(:,j)=R(1:2,1:2)*pos(:,j)+obj.start_pos(1:2);
            end

            plot(pos(1,:),pos(2,:),'-o','MarkerSize',2)
            
            hold on
            axis equal
            
            g_exp(1:3,1:3)=rotz(obj.start_pos(3))*g_exp(1:3,1:3);
            g_exp(1:2,4)=pos(:,end);

            new_xy=g_exp(1:2,1:2)*[1 0;0 1];
            new_x=new_xy(1,:);
            new_y=new_xy(2,:);
            End_pos=g_exp(1:2,4);
            location_x=[End_pos(1),End_pos(1)];
            location_y=[End_pos(2),End_pos(2)];
            quiver(location_x,location_y,new_x,new_y,0.05*obj.Ltotal)
       end
       

       %% ������е������
       function cal_pos_all(obj)
            theta_solve=obj.conv_theta;
            Theta=zeros(1,length(theta_solve));
            for k=1:length(Theta)
                Theta(k)=theta_solve(k);
            end

            g_exp=obj.conv_g;

            %�õ������˵�λ��
            pos=[0;0];
            pos=[pos,[obj.conv_q_all(1,1);obj.conv_q_all(2,1)]];
            for i =2:length(Theta)
                g0i=[eye(3),[(i-0.5)*obj.seg_length;0;0];0,0,0,1];
                gsti=obj.exp_fkine(obj.conv_w_all(:,1:i-1),obj.conv_q_all(:,1:i-1),g0i,Theta(:,1:i-1));
                pos=[pos,gsti(1:2,4)];
            end
            pos=[pos,g_exp(1:2,4)];
            
            obj.conv_pos_all=pos;
            
            for j=1:length(pos)
                R=rotz(obj.start_pos(3));
                pos(:,j)=R(1:2,1:2)*pos(:,j)+obj.start_pos(1:2);
            end
            
            obj.pos_all=pos;
                 
       end
        
    end
end