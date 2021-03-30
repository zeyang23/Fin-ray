classdef finray_contact < handle
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
        % ��ߵĲ���
        wid_A;
        thi_A;
        E_A;
        % �ұߵĲ���
        wid_B;
        thi_B;
        E_B;
        
        % �ڲ�����Լ��
        constraint_number;
        
        % 3*m���� mΪ����Լ���ĸ���
        constraint_ratio;  % Լ������+��߱���+�ұ߱��� 
        constraint_index; % Լ������+��߶�Ӧ�ؽڱ��+�ұ�
        
        
        % ��߽Ӵ�
        % �Ӵ�������
        A_contact_index;
        A_contact_force;
        
        
        % ���㵼�������planar_nR���󣬶�������
        RodA;
        RodB;
        
        Rod_A_contact;
        
        A_constraint_array;
        B_constraint_array;
        
        
        % Բ��λ��
        xc;
        yc;
        rc;
        
    end
    
    methods
        function obj = finray_contact(finray_info)
            % finray_info��һ���ṹ�壬�����˱�Ҫ����Ϣ
            
            % ��finray_info��ȡ��Ϣ
            obj.xc=finray_info.xc;
            obj.yc=finray_info.yc;
            
            obj.rc=finray_info.rc;
            
            obj.pA=finray_info.pA;
            obj.pB=finray_info.pB;
            
            obj.LA=finray_info.LA;
            obj.LB=finray_info.LB;
            
            obj.nA=finray_info.nA;
            obj.nB=finray_info.nB;
            
            obj.psi=finray_info.psi;
            
            obj.wid_A=finray_info.wid_A;
            obj.wid_B=finray_info.wid_B;
            
            obj.thi_A=finray_info.thi_A;
            obj.thi_B=finray_info.thi_B;
            
            obj.E_A=finray_info.E_A;
            obj.E_B=finray_info.E_B;
            
            obj.constraint_ratio=finray_info.constraint_ratio;
            
            obj.A_contact_index=finray_info.A_contact_index;
            obj.A_contact_force=finray_info.A_contact_force;
            
            
            % ��ʼ������Ϣ
            obj.RodA=planar_nR(obj.E_A,obj.LA,obj.wid_A,obj.thi_A,obj.nA,[0;0;0]);
            obj.RodB=planar_nR(obj.E_B,obj.LB,obj.wid_B,obj.thi_B,obj.nB,[0;0;0]);
            
            
            obj.A_constraint_array=planar_nR.empty;
            obj.B_constraint_array=planar_nR.empty;
            
            obj.Rod_A_contact=planar_nR(obj.E_A,obj.A_contact_index/obj.nA*obj.LA,obj.wid_A,obj.thi_A,obj.A_contact_index,[0;0;0]);
            
            
            
            obj.constraint_number=size(obj.constraint_ratio,1);
            obj.constraint_index=zeros(size(obj.constraint_ratio));
            obj.constraint_index(:,1)=obj.constraint_ratio(:,1);

            for i=1:obj.constraint_number
                kai=fix(obj.constraint_ratio(i,2)*obj.nA);
                obj.constraint_index(i,2)=kai;
                obj.A_constraint_array(i)=planar_nR(obj.E_A,kai/obj.nA*obj.LA,obj.wid_A,obj.thi_A,kai,[0;0;0]);

                kbi=fix(obj.constraint_ratio(i,3)*obj.nB);
                obj.constraint_index(i,3)=kbi;
                obj.B_constraint_array(i)=planar_nR(obj.E_B,kbi/obj.nB*obj.LA,obj.wid_B,obj.thi_B,kbi,[0;0;0]);
            end
        end

        
        
        function [r,J]=cal_balance(obj,x)
            alpha=obj.pA(3);
            xA=obj.pA(1);
            yA=obj.pA(2);
            
            beta=obj.pB(3);
            xB=obj.pB(1);
            yB=obj.pB(2);
            
            alpha_degree=alpha/pi*180;
            beta_degree=beta/pi*180;

            A=rotz(-alpha_degree)*rotz(beta_degree);
            b=rotz(-alpha_degree)*[xB-xA;yB-yA;beta-alpha-obj.psi];
            B=-rotz(beta_degree-alpha_degree);
            
            
            
            % ��x��ȡ����
            
            thetaA=x(1:obj.nA);
            thetaB=x(obj.nA+1:obj.nA+obj.nB);
            FB=x(obj.nA+obj.nB+1:obj.nA+obj.nB+3);

            FA=B*FB;

            obj.RodA.theta=thetaA;
            obj.RodA.F=FA;
            obj.RodB.theta=thetaB;
            obj.RodB.F=FB;

            obj.RodA.update;
            obj.RodB.update;
            
            
            obj.rc=x(end);
            
            theta_A_contact=thetaA(1:obj.A_contact_index);
            obj.A_contact_force=x(end-1);
            
            obj.Rod_A_contact.theta=theta_A_contact;
            obj.Rod_A_contact.F=zeros(3,1);
            obj.Rod_A_contact.update;
            
            p_contact=obj.Rod_A_contact.pe;
            obj.Rod_A_contact.F=[obj.A_contact_force*sin(p_contact(3));-obj.A_contact_force*cos(p_contact(3));0];
            obj.Rod_A_contact.update;
            
            
            % �и���Լ�������    

            r=zeros(obj.nA+obj.nB+5+2*obj.constraint_number,1);
            J=zeros(obj.nA+obj.nB+5+2*obj.constraint_number);

            r1=obj.RodA.pe-A*obj.RodB.pe-b;
            r2=obj.RodA.K_theta*thetaA-transpose(obj.RodA.Jacobian)*FA;
            r3=obj.RodB.K_theta*thetaB-transpose(obj.RodB.Jacobian)*FB;

            r(1:obj.nA+obj.nB+3)=[r1;r2;r3];

            J(1:3,1:obj.nA)=obj.RodA.Jacobian;
            J(1:3,obj.nA+1:obj.nA+obj.nB)=-A*obj.RodB.Jacobian;

            J(4:obj.nA+3,1:obj.nA)=obj.RodA.K_theta-obj.RodA.partial;
            J(4:obj.nA+3,obj.nA+obj.nB+1:obj.nA+obj.nB+3)=-transpose(obj.RodA.Jacobian)*B;

            J(obj.nA+4:obj.nA+obj.nB+3,obj.nA+1:obj.nA+obj.nB)=obj.RodB.K_theta-obj.RodB.partial;
            J(obj.nA+4:obj.nA+obj.nB+3,obj.nA+obj.nB+1:obj.nA+obj.nB+3)=-transpose(obj.RodB.Jacobian);
            
            
            % �е�λ�÷���
            p_contact=obj.Rod_A_contact.pe;
            r4=[obj.rc*sin(p_contact(3)+obj.pA(3))-(p_contact(1)*cos(obj.pA(3))-p_contact(2)*sin(obj.pA(3)))+obj.xc-obj.pA(1);
                -obj.rc*cos(p_contact(3)+obj.pA(3))-(p_contact(1)*sin(obj.pA(3))+p_contact(2)*cos(obj.pA(3)))+obj.yc-obj.pA(2)];
            
            r(end-1:end)=r4;
            

            % �������Լ������
            for i=1:obj.constraint_number
                Lcon=obj.constraint_index(i,1);

                fcon=x(obj.nA+obj.nB+3+2*i-1);
                gamma=x(obj.nA+obj.nB+3+2*i);

                obj.A_constraint_array(i).theta=thetaA(1:obj.constraint_index(i,2));
                obj.A_constraint_array(i).F=-rotz(-alpha_degree)*[fcon*cos(gamma);fcon*sin(gamma);0];
                obj.A_constraint_array(i).update;

                obj.B_constraint_array(i).theta=thetaB(1:obj.constraint_index(i,3));
                obj.B_constraint_array(i).F=rotz(-beta_degree)*[fcon*cos(gamma);fcon*sin(gamma);0];
                obj.B_constraint_array(i).update;

                lambdaAi=zeros(obj.constraint_index(i,2),obj.nA);
                lambdaAi(1:obj.constraint_index(i,2),1:obj.constraint_index(i,2))=diag(ones(obj.constraint_index(i,2),1));

                lambdaBi=zeros(obj.constraint_index(i,3),obj.nB);
                lambdaBi(1:obj.constraint_index(i,3),1:obj.constraint_index(i,3))=diag(ones(obj.constraint_index(i,3),1));

                % ����������������
                r(4:obj.nA+3)=r(4:obj.nA+3)-transpose(lambdaAi)*transpose(obj.A_constraint_array(i).Jacobian)*obj.A_constraint_array(i).F;
                r(obj.nA+4:obj.nA+obj.nB+3)=r(obj.nA+4:obj.nA+obj.nB+3)-transpose(lambdaBi)*transpose(obj.B_constraint_array(i).Jacobian)*obj.B_constraint_array(i).F;

                % �����ڲ�Լ������
                pkai=obj.A_constraint_array(i).pe;
                pkbi=obj.B_constraint_array(i).pe;
                r(obj.nA+obj.nB+3+2*i-1:obj.nA+obj.nB+3+2*i)=[pkai(1);pkai(2)]-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*[pkbi(1);pkbi(2)]+...
                [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[xA-xB+Lcon*cos(gamma);yA-yB+Lcon*sin(gamma)];

                % ����A������ƽ�ⷽ�̵ĵ�����3��
                J(4:obj.nA+3,1:obj.nA)=J(4:obj.nA+3,1:obj.nA)-transpose(lambdaAi)*obj.A_constraint_array(i).partial*lambdaAi;
                J(4:obj.nA+3,obj.nA+obj.nB+3+2*i-1)=-transpose(lambdaAi)*transpose(obj.A_constraint_array(i).Jacobian)*(-rotz(-alpha_degree)*[cos(gamma);sin(gamma);0]);
                J(4:obj.nA+3,obj.nA+obj.nB+3+2*i)=-transpose(lambdaAi)*transpose(obj.A_constraint_array(i).Jacobian)*(-rotz(-alpha_degree)*[-fcon*sin(gamma);fcon*cos(gamma);0]);

                % ����B������ƽ�ⷽ�̵ĵ�����3��
                J(obj.nA+4:obj.nA+obj.nB+3,obj.nA+1:obj.nA+obj.nB)=J(obj.nA+4:obj.nA+obj.nB+3,obj.nA+1:obj.nA+obj.nB)-transpose(lambdaBi)*obj.B_constraint_array(i).partial*lambdaBi;
                J(obj.nA+4:obj.nA+obj.nB+3,obj.nA+obj.nB+3+2*i-1)=-transpose(lambdaBi)*transpose(obj.B_constraint_array(i).Jacobian)*(rotz(-beta_degree)*[cos(gamma);sin(gamma);0]);
                J(obj.nA+4:obj.nA+obj.nB+3,obj.nA+obj.nB+3+2*i)=-transpose(lambdaBi)*transpose(obj.B_constraint_array(i).Jacobian)*(rotz(-beta_degree)*[-fcon*sin(gamma);fcon*cos(gamma);0]);

                % �����ڲ�ƽ�ⷽ�̵�3���
                tempA=obj.A_constraint_array(i).Jacobian*lambdaAi;
                J(obj.nA+obj.nB+3+2*i-1:obj.nA+obj.nB+3+2*i,1:obj.nA)=tempA(1:2,:);
                tempB=obj.B_constraint_array(i).Jacobian*lambdaBi;
                J(obj.nA+obj.nB+3+2*i-1:obj.nA+obj.nB+3+2*i,obj.nA+1:obj.nA+obj.nB)=-[cos(beta-alpha) -sin(beta-alpha); sin(beta-alpha) cos(beta-alpha)]*tempB(1:2,:);
                J(obj.nA+obj.nB+3+2*i-1:obj.nA+obj.nB+3+2*i,obj.nA+obj.nB+3+2*i)=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*[-Lcon*sin(gamma);Lcon*cos(gamma)];

            end


            % ����Ӵ�������
            lambdaF=zeros(obj.A_contact_index,obj.nA);
            lambdaF(1:obj.A_contact_index,1:obj.A_contact_index)=diag(ones(obj.A_contact_index,1));
            r(4:obj.nA+3)=r(4:obj.nA+3)-transpose(lambdaF)*transpose(obj.Rod_A_contact.Jacobian)*obj.Rod_A_contact.F;

            J(4:obj.nA+3,1:obj.nA)=J(4:obj.nA+3,1:obj.nA)...
                -transpose(lambdaF)*obj.Rod_A_contact.partial*lambdaF...
                -transpose(lambdaF)*transpose(obj.Rod_A_contact.Jacobian)*[obj.A_contact_force*cos(p_contact(3));obj.A_contact_force*sin(p_contact(3));0]*ones(1,obj.A_contact_index)*lambdaF;

            J(4:obj.nA+3,end-1)=-transpose(lambdaF)*transpose(obj.Rod_A_contact.Jacobian)*[sin(p_contact(3));-cos(p_contact(3));0];
            
            
            % �е�λ�÷��̵ĵ���
            % ��theta_a�ĵ���
            temp=obj.Rod_A_contact.Jacobian*lambdaF;
            J(end-1:end,1:obj.nA)=[obj.rc*cos(p_contact(3)+obj.pA(3))*[ones(1,obj.A_contact_index),zeros(1,obj.nA-obj.A_contact_index)];
                                   obj.rc*sin(p_contact(3)+obj.pA(3))*[ones(1,obj.A_contact_index),zeros(1,obj.nA-obj.A_contact_index)]]-[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]*temp(1:2,:);
            
            % �԰뾶�ĵ���
            J(end-1:end,end)=[sin(p_contact(3)+obj.pA(3));-cos(p_contact(3)+obj.pA(3))];
            
        end
        
        
        function plot_state(obj,x)
            
            % ����A��B����״
            thetaA=x(1:obj.nA);
            thetaB=x(obj.nA+1:obj.nA+obj.nB);
            
            obj.RodA.theta=thetaA;
            obj.RodA.cal_pe;
            obj.RodA.cal_posall;
            plot_abs_pos(obj.RodA.pos_all,obj.pA(3),[obj.pA(1),obj.pA(2)]);
            hold on
            
            obj.RodB.theta=thetaB;
            obj.RodB.cal_pe;
            obj.RodB.cal_posall;
            plot_abs_pos(obj.RodB.pos_all,obj.pB(3),[obj.pB(1),obj.pB(2)]);
            
            
            for i=1:obj.constraint_number
                obj.A_constraint_array(i).theta=thetaA(1:obj.constraint_index(i,2));
                obj.A_constraint_array(i).cal_pe;

                obj.B_constraint_array(i).theta=thetaB(1:obj.constraint_index(i,3));
                obj.B_constraint_array(i).cal_pe;

                pka=obj.A_constraint_array(i).pe;
                PA(1)=pka(1)*cos(obj.pA(3))-pka(2)*sin(obj.pA(3))+obj.pA(1);
                PA(2)=pka(1)*sin(obj.pA(3))+pka(2)*cos(obj.pA(3))+obj.pA(2);

                pkb=obj.B_constraint_array(i).pe;
                PB(1)=pkb(1)*cos(obj.pB(3))-pkb(2)*sin(obj.pB(3))+obj.pB(1);
                PB(2)=pkb(1)*sin(obj.pB(3))+pkb(2)*cos(obj.pB(3))+obj.pB(2);

                plot([PA(1) PB(1)],[PA(2) PB(2)])
            end

            
            scale=0.05;
            % ����A������λ��
            
        end
    end
end