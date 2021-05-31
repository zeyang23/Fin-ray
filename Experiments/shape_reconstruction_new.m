function [X,A_abs_pos,B_abs_pos,sensor1_abs_pos,sensor2_abs_pos,contact_abs_pos,rigid_pos_left,rigid_pos_right]=shape_reconstruction_new(X0,U1,U2)

    % ���ݴ�����ʾ��������н�

    coff_1=[1.186256198447459,0.056997358385406];
    coff_2=[1.151462241145653,0.027617725599564];


    U1_0=3.15-1.46;
    U2_0=1.46;

    U1_delta=U1-U1_0;
    U2_delta=U2-U2_0;


    DeltaA=polyval(coff_1,U1_delta);
    DeltaB=polyval(coff_2,U2_delta);


    lambda1=0.45;
    lambda2=0.45;



    xA=0;
    yA=0;

    xB=73.24e-3;
    yB=0;


    alpha_degree=90;
    beta_degree=116;
    psi_degree=beta_degree-alpha_degree;

    psi=psi_degree/180*pi;
    alpha=alpha_degree/180*pi;
    beta=beta_degree/180*pi;


    nA=50;
    nB=50;


    LA=(xB-xA)/tan(psi);
    LB=(xB-xA)/sin(psi);

    wid=28e-3;
    thi=0.2e-3;
    E=213e9;


    constraint_ratio=[];

    La_1=27e-3;
    La_2=54e-3;
    La_3=81e-3;
    La_4=108e-3;

    Lb_1=45.67e-3;
    Lb_2=78.69e-3;
    Lb_3=107.5e-3;
    Lb_4=132.27e-3;

    Lcon1=sqrt((LA-La_1)^2+(LB-Lb_1)^2-2*(LA-La_1)*(LB-Lb_1)*cos(psi));
    Lcon2=sqrt((LA-La_2)^2+(LB-Lb_2)^2-2*(LA-La_2)*(LB-Lb_2)*cos(psi));
    Lcon3=sqrt((LA-La_3)^2+(LB-Lb_3)^2-2*(LA-La_3)*(LB-Lb_3)*cos(psi));
    Lcon4=sqrt((LA-La_4)^2+(LB-Lb_4)^2-2*(LA-La_4)*(LB-Lb_4)*cos(psi));

    constraint_ratio=[Lcon1,La_1/LA,Lb_1/LB;
                      Lcon2,La_2/LA,Lb_2/LB;
                      Lcon3,La_3/LA,Lb_3/LB;
                      Lcon4,La_4/LA,Lb_4/LB];

    A_force_ratio=[];
    A_force_ratio=[0,lambda1;
                   0,lambda2];


    B_force_ratio=[];


    % ���������趨�Ĳ������ɲ����ṹ��
    finray_info=struct();

    finray_info.pA=[xA;yA;alpha];
    finray_info.pB=[xB;yB;beta];

    finray_info.LA=LA;
    finray_info.LB=LB;

    finray_info.nA=nA;
    finray_info.nB=nB;

    finray_info.psi=psi;

    finray_info.wid_A=wid;
    finray_info.wid_B=wid;

    finray_info.thi_A=thi;
    finray_info.thi_B=thi;

    finray_info.E_A=LA/135e-3*E;
    finray_info.E_B=LB/156e-3*E;

    finray_info.constraint_ratio=constraint_ratio;

    finray_info.A_force_ratio=A_force_ratio;

    finray_info.B_force_ratio=B_force_ratio;


    Finray_solve=finray_force(finray_info);


    % �������ĳ���
    sensor_length=36e-3;

    % ������������λ��
    sensor1_center=133.9e-3;
    sensor2_center=93.1e-3;


    p1=fix((sensor1_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
    q1=fix((sensor1_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;

    p2=fix((sensor2_center-sensor_length/2-LB/nB/2)/(LB/nB))+1;
    q2=fix((sensor2_center+sensor_length/2-LB/nB/2)/(LB/nB))+1;




    x0=zeros(nA+nB+3+2*Finray_solve.constraint_number+2,1);

    f=@(x) cal_all_balance(x,Finray_solve,p1,q1,p2,q2,DeltaA,DeltaB);

    options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'CheckGradient',false,'Display','off','Algorithm','levenberg-marquardt');

    [X,fval,exitflag,output] = fsolve(f,X0,options);

    
    figure
    [A_abs_pos,B_abs_pos]=Finray_solve.plot_state_sr(X(1:end-2));



    Finray_solve.RodB.cal_posall;
    sensor1_pos=Finray_solve.RodB.pos_all(1+p1:1+q1,:);
    sensor2_pos=Finray_solve.RodB.pos_all(1+p2:1+q2,:);

    
    sensor1_abs_pos=plot_sensor_pos(sensor1_pos,beta,[xB,yB]);
    sensor2_abs_pos=plot_sensor_pos(sensor2_pos,beta,[xB,yB]);
    
    
    
    contact_center_length=0.45*LA;
    contact_length=0.3*LA;
    
    p_contact=fix((contact_center_length-contact_length/2-LA/nA/2)/(LA/nA))+1;
    q_contact=fix((contact_center_length+contact_length/2-LA/nA/2)/(LA/nA))+1;
    
    
    Finray_solve.RodA.cal_posall;
    contact_pos=Finray_solve.RodA.pos_all(1+p_contact:1+q_contact,:);
    
    contact_abs_pos=plot_sensor_pos(contact_pos,alpha,[xA,yA]);
    
    
    rigid_pos_left=[];
    rigid_pos_right=[];
    
    for i=1:Finray_solve.constraint_number
        Finray_solve.A_constraint_array(i).theta=Finray_solve.RodA.theta(1:Finray_solve.constraint_index(i,2));
        Finray_solve.A_constraint_array(i).cal_pe;

        Finray_solve.B_constraint_array(i).theta=Finray_solve.RodB.theta(1:Finray_solve.constraint_index(i,3));
        Finray_solve.B_constraint_array(i).cal_pe;

        pka=Finray_solve.A_constraint_array(i).pe;
        PA(1)=pka(1)*cos(Finray_solve.pA(3))-pka(2)*sin(Finray_solve.pA(3))+Finray_solve.pA(1);
        PA(2)=pka(1)*sin(Finray_solve.pA(3))+pka(2)*cos(Finray_solve.pA(3))+Finray_solve.pA(2);

        pkb=Finray_solve.B_constraint_array(i).pe;
        PB(1)=pkb(1)*cos(Finray_solve.pB(3))-pkb(2)*sin(Finray_solve.pB(3))+Finray_solve.pB(1);
        PB(2)=pkb(1)*sin(Finray_solve.pB(3))+pkb(2)*cos(Finray_solve.pB(3))+Finray_solve.pB(2);

        rigid_pos_left(i,:)=[PA(1) PA(2)];
        rigid_pos_right(i,:)=[PB(1) PB(2)];

    end
    
end


function [r,J]=cal_all_balance(x,Finray,p1,q1,p2,q2,Delta1,Delta2)

    Finray.A_force_index(1,1)=x(end-1);
    Finray.A_force_index(2,1)=x(end);
    
    Finray.A_force_ratio(1,1)=x(end-1);
    Finray.A_force_ratio(2,1)=x(end);
    
    [r1,J1]=Finray.cal_balance(x(1:end-2));
    
    r=zeros(length(r1)+2,1);
    J=zeros(length(r1)+2);
    
    r(1:length(r1))=r1;
    J(1:length(r1),1:length(r1))=J1;
    
    Delta1_now=sum(Finray.RodB.theta(p1:q1));
    Delta2_now=sum(Finray.RodB.theta(p2:q2));
    
    r(length(r1)+1)=Delta1_now-Delta1;
    r(length(r1)+2)=Delta2_now-Delta2;
    
    D1=zeros(1,Finray.RodB.n_seg);
    for i=p1:q1
        D1(i)=1;
    end
    
    D2=zeros(1,Finray.RodB.n_seg);
    for i=p2:q2
        D2(i)=1;
    end
    
    
    
    JF1=Finray.A_force_array(1).Jacobian;
    JF2=Finray.A_force_array(2).Jacobian;
    
    lambdaF1=zeros(Finray.A_force_index(1,2),Finray.nA);
    lambdaF1(1:Finray.A_force_index(1,2),1:Finray.A_force_index(1,2))=diag(ones(Finray.A_force_index(1,2),1));
    
    lambdaF2=zeros(Finray.A_force_index(2,2),Finray.nA);
    lambdaF2(1:Finray.A_force_index(2,2),1:Finray.A_force_index(2,2))=diag(ones(Finray.A_force_index(2,2),1));
    
    pk_F1=Finray.A_force_array(1).pe;
    pk_F2=Finray.A_force_array(2).pe;
    
    J(4:Finray.RodA.n_seg+3,end-1)=-transpose(lambdaF1)*transpose(JF1)*[sin(pk_F1(3));-cos(pk_F1(3));0];
    J(4:Finray.RodA.n_seg+3,end)=-transpose(lambdaF2)*transpose(JF2)*[sin(pk_F2(3));-cos(pk_F2(3));0];
    
    J(length(r1)+1,Finray.nA+1:Finray.nA+Finray.nB)=D1;
    J(length(r1)+2,Finray.nA+1:Finray.nA+Finray.nB)=D2;
end


function abs_pos=plot_sensor_pos(pos_all,alpha,pA)
    xA=pA(1);
    yA=pA(2);
    abs_pos=zeros(size(pos_all));
    n=size(pos_all,1);
    for i=1:n
        abs_pos(i,1)=pos_all(i,1)*cos(alpha)-pos_all(i,2)*sin(alpha)+xA;
        abs_pos(i,2)=pos_all(i,1)*sin(alpha)+pos_all(i,2)*cos(alpha)+yA;
    end
%     plot(abs_pos(:,1),abs_pos(:,2),'o','MarkerSize',4,'LineWidth',2);
%     axis equal
end