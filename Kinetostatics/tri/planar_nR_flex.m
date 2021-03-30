% 不使用指数坐标，仅使用三角函数进行计算

% 长度可变的柔性板，用于处理相切情况

classdef planar_nR_flex < handle
    
    properties   
        Ltotal;
        n_seg;
        seg_length;
        
        wid; % 柔性板的截面长度
        thi; % 柔性板的截面宽度
        cross_section; % 柔性板的横截面积
        Iz; % 截面惯性矩        
        E; % 杨氏模量
        
        K_theta;
        
        F % 3*1
        
        theta; % n*1
        
        pe;
        Jacobian;
        
        Hessian_xe;
        Hessian_ye;
        
        partial; % 计算那项奇怪的偏导数
        
        pe_desired % 3*1
        
        pos_all % (n+2)*2
    end
    
    methods
        
        function item = planar_nR_flex(E,l,w,h,n,pdes)
           item.E=E;
           item.Ltotal=l;
           item.wid=w;
           item.thi=h;
           item.n_seg=n;
           item.pe_desired=pdes;
           
           item.cross_section=item.wid*item.thi;
           item.Iz=1/12*item.thi.^3*item.wid;
           item.seg_length=item.Ltotal/item.n_seg;
           item.K_theta=diag((item.E*item.Iz/item.seg_length)*ones(n,1));
           
           item.theta=zeros(item.n_seg,1);
           item.F=zeros(3,1);
           item.update;
        end
        
        function change_length(obj,L_new)
            obj.Ltotal=L_new;
            obj.seg_length=obj.Ltotal/obj.n_seg;
            obj.K_theta=diag((obj.E*obj.Iz/obj.seg_length)*ones(obj.n_seg,1));
            obj.update;
        end
        
        function cal_pe(obj)
            xe=0;
            ye=0;
            phie=sum(obj.theta);
            
            for i=1:obj.n_seg-1
                xe=xe+obj.seg_length*cos(sum(obj.theta(1:i)));
                ye=ye+obj.seg_length*sin(sum(obj.theta(1:i)));
            end
            
            xe=xe+obj.seg_length/2+obj.seg_length*cos(sum(obj.theta))/2;
            ye=ye+obj.seg_length*sin(sum(obj.theta))/2;
            obj.pe=[xe;ye;phie];
        end
        
        function pk=cal_pk(obj,k)
            xk=0;
            yk=0;
            phik=sum(obj.theta(1:k));
            
            for i=1:k
                xk=xk+obj.seg_length*cos(sum(obj.theta(1:i)));
                yk=yk+obj.seg_length*sin(sum(obj.theta(1:i)));
            end
            
            xk=xk+obj.seg_length/2;
            pk=[xk;yk;phik];
        end
        
        function cal_Jacobian(obj)
            J=zeros(3,obj.n_seg);
            for i=1:obj.n_seg-1
                J(:,i)=zeros(3,1);
                J(3,i)=1;
                for j=i:obj.n_seg-1
                    J(1,i)=J(1,i)+(-1)*obj.seg_length*sin(sum(obj.theta(1:j)));
                    J(2,i)=J(2,i)+obj.seg_length*cos(sum(obj.theta(1:j)));
                end
                J(1,i)=J(1,i)+(-1)*obj.seg_length*sin(sum(obj.theta))/2;
                J(2,i)=J(2,i)+obj.seg_length*cos(sum(obj.theta))/2;
            end
            
            J(1,end)=(-1)*obj.seg_length*sin(sum(obj.theta))/2;
            J(2,end)=obj.seg_length*cos(sum(obj.theta))/2;
            J(3,end)=1;
            obj.Jacobian=J;
        end
        
        function cal_Hessian(obj)
            Hx=zeros(obj.n_seg);
            for i=1:obj.n_seg-1
                temp=0;
                for j=i:obj.n_seg-1
                    temp=temp+(-1)*obj.seg_length*cos(sum(obj.theta(1:j)));
                end
                temp=temp-obj.seg_length*cos(sum(obj.theta))/2;
                Hx(1:i,i)=temp*ones(i,1);
                
                for j=i+1:obj.n_seg-1
                    for k=j:obj.n_seg-1
                        Hx(j,i)=Hx(j,i)+(-1)*obj.seg_length*cos(sum(obj.theta(1:k)));
                    end
                end
                Hx(i+1:end,i)=Hx(i+1:end,i)+(-obj.seg_length*cos(sum(obj.theta))/2)*ones(obj.n_seg-i,1);
            end
            
            Hx(:,end)=-1/2*obj.seg_length*cos(sum(obj.theta))*ones(obj.n_seg,1);
            obj.Hessian_xe=Hx;
            
            Hy=zeros(obj.n_seg);
            for i=1:obj.n_seg-1
                temp=0;
                for j=i:obj.n_seg-1
                    temp=temp+(-1)*obj.seg_length*sin(sum(obj.theta(1:j)));
                end
                temp=temp-obj.seg_length*sin(sum(obj.theta))/2;
                Hy(1:i,i)=temp*ones(i,1);
                
                for j=i+1:obj.n_seg-1
                    for k=j:obj.n_seg-1
                        Hy(j,i)=Hy(j,i)+(-1)*obj.seg_length*sin(sum(obj.theta(1:k)));
                    end
                end
                Hy(i+1:end,i)=Hy(i+1:end,i)+(-obj.seg_length*sin(sum(obj.theta))/2)*ones(obj.n_seg-i,1);
            end
            
            Hy(:,end)=-1/2*obj.seg_length*sin(sum(obj.theta))*ones(obj.n_seg,1);
            obj.Hessian_ye=Hy;
        end
        
        function cal_partial(obj)
            % 注意，调用这个函数前记得先更新Hessian矩阵
            
            temp=zeros(obj.n_seg);
            for i=1:obj.n_seg
                for j=1:obj.n_seg
                    temp(i,j)=transpose(obj.F)*[obj.Hessian_xe(j,i);obj.Hessian_ye(j,i);0];
                end
            end
            
            obj.partial=temp;
        end
            
        function update(obj)
            obj.cal_pe;
            obj.cal_Jacobian;
            obj.cal_Hessian;
            obj.cal_partial;
        end
        
        function r=cal_r(obj)
            r=zeros(obj.n_seg+3,1);
            r(1:3)=obj.pe-obj.pe_desired;
            r(4:end)=obj.K_theta*obj.theta-transpose(obj.Jacobian)*obj.F;
        end
        
        function J=cal_rdot(obj)
            J=zeros(obj.n_seg+3);
            J(1:3,1:obj.n_seg)=obj.Jacobian;
            J(4:obj.n_seg+3,end-2:end)=-transpose(obj.Jacobian);
            
            temp=zeros(obj.n_seg);
            for i=1:obj.n_seg
                for j=1:obj.n_seg
                    temp(i,j)=transpose(obj.F)*[obj.Hessian_xe(j,i);obj.Hessian_ye(j,i);0];
                end
            end
            
            J(4:end,1:obj.n_seg)=obj.K_theta-0*temp;
        end
        
        function Newton(obj,TOL)
            % 牛顿法求解器
            k=1;
            while(1)
                if k>500
                    error("can not converge")
                end
                J=obj.cal_rdot;
                b=obj.cal_r;
                delta=-pinv(J)*b;
                
                obj.theta=obj.theta+delta(1:obj.n_seg);
                obj.F=obj.F+delta(end-2:end);
                
                obj.update;
                
                k=k+1;
                if(norm(obj.cal_r)<TOL)
                    fprintf('Newton Method converge: iteration = %d\n',k-1)
                    fprintf('norm(e) = %E\n',norm(obj.cal_r))
                    break;
                end
            end
        end
        
        function cal_posall(obj)
            posall=zeros(obj.n_seg+2,2);
            posall(2,:)=[1/2*obj.seg_length,0];
            posall(end,:)=[obj.pe(1),obj.pe(2)];
            
            for i=1:obj.n_seg-1
                posall(i+2,1)=posall(2,1);
                posall(i+2,2)=posall(2,2);
                for j=1:i
                    posall(i+2,1)=posall(i+2,1)+obj.seg_length*cos(sum(obj.theta(1:j)));
                    posall(i+2,2)=posall(i+2,2)+obj.seg_length*sin(sum(obj.theta(1:j)));
                end
            end
            
            obj.pos_all=posall;
        end
        
        function plot_all(obj)
            obj.cal_posall;
            pos=obj.pos_all;
            plot(pos(:,1),pos(:,2),'-o','MarkerSize',2);
            axis equal
        end
    end
end