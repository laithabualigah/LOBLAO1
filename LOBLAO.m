

% enahnced Aquila Optimizer
function [Best_FF,Best_P,Conv_curve]=LOBLAO(N,M_Iter,LB,UB,Dim,k,D)
display('LOBLAO Working');
%Two variables to keep the positions and the fitness value of the best-obtained solution

% Best_P=zeros(1,Dim);
Best_FF=inf;
Conv_curve=zeros(1,M_Iter);

%Initialize the positions of solution
X=initialization(N,Dim,UB,LB);
Xnew=X;
% Ffun=zeros(1,size(X,1));% (fitness values)
% Ffun_new=zeros(1,size(Xnew,1));% (fitness values)

MOP_Max=1;
MOP_Min=0.2;
C_Iter=1;
Alpha=5;
Mu=0.499;


for i=1:size(X,1)
    Ffun(1,i)=fitn(k,X(i,:),D);
    if Ffun(1,i)<Best_FF
        Best_FF=Ffun(1,i);
        Best_P=X(i,:);
    end
end
    
   


alpha=0.1;
delta=0.1;

while C_Iter<M_Iter+1  %Main loop
    for i=1:size(X,1)
        F_UB=X(i,:)>UB;
        F_LB=X(i,:)<LB;
        X(i,:)=(X(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
        Ffun_new(1,i)=fitn(k,Xnew(i,:),D);
        if Ffun_new(1,i)<Best_FF
            Best_FF=fitn(1,i);
            Best_P=X(i,:);
        end
    end

    
    G2=2*rand()-1; % Eq. (16)
    G1=2*(1-(C_Iter/M_Iter));  % Eq. (17)
    to = 1:Dim;
    u = .0265;
    r0 = 10;
    r = r0 +u*to;
    omega = .005;
    phi0 = 3*pi/2;
    phi = -omega*to+phi0;
    x = r .* sin(phi);  % Eq. (9)
    y = r .* cos(phi); % Eq. (10)
    QF=C_Iter^((2*rand()-1)/(1-M_Iter)^2); % Eq. (15)
        %-------------------------------------------------------------------------------------
    for i=1:size(X,1)
        %-------------------------------------------------------------------------------------
        if C_Iter<=(2/3)*M_Iter
            if rand <0.5
                Xnew(i,:)=Best_P(1,:)*(1-C_Iter/M_Iter)+(mean(X(i,:))-Best_P(1,:))*rand(); % Eq. (3) and Eq. (4)
         Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=round((Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB);
                Ffun_new(1,i)=fitn(k,Xnew(i,:),D);
                if Ffun_new(1,i)<Best_FF
                    X(i,:)=Xnew(i,:);
                    Best_FF=Ffun_new(1,i);
                end
            else
                %-------------------------------------------------------------------------------------
                Xnew(i,:)=Best_P(1,:).*Levy1(Dim)+X((floor(N*rand()+1)),:)+(y-x)*rand;       % Eq. (5)
                        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=round((Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB);
                Ffun_new(1,i)=fitn(k,Xnew(i,:),D);
                if Ffun_new(1,i)<Best_FF
                    X(i,:)=Xnew(i,:);
                    Best_FF=Ffun_new(1,i);
                end
            end
            %-------------------------------------------------------------------------------------
        else
            if rand<0.5
                Xnew(i,:)=(Best_P(1,:)-mean(X))*alpha-rand+((UB-LB)*rand+LB)*delta;   % Eq. (13)
                        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=round((Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB);
                Ffun_new(1,i)=fitn(k,Xnew(i,:),D);
                if Ffun_new(1,i)<Best_FF
                    X(i,:)=Xnew(i,:);
                    Best_FF=Ffun_new(1,i);
                end
            else
                %-------------------------------------------------------------------------------------
                Xnew(i,:)=QF*Best_P(1,:)-(G2*X(i,:)*rand)-G1.*Levy1(Dim)+rand*G2; % Eq. (14)
                        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=round((Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB);
              Ffun_new(1,i)=fitn(k,Xnew(i,:),D);
                if Ffun_new(1,i)<Best_FF
                    X(i,:)=Xnew(i,:);
                    Best_FF=Ffun_new(1,i);
                end
            end
            
        end
    end
 

    
    
    
                % Mutation Search Strategy (MSS)
            if rand < 1
                mutate_index = randi(Dim);
                Xnew(i, mutate_index) = LB + rand*(UB-LB);
            end
                                    Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=round((Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB);
              Ffun_new(1,i)=fitn(k,Xnew(i,:),D);
                if Ffun_new(1,i)<Best_FF
                    X(i,:)=Xnew(i,:);
                    Best_FF=Ffun_new(1,i);
                end
            
            

                
             for i = 1:N
            % Opposition-Based Learning (OBL)
            X_OBL = LB + UB - X(i, :);
            Ffun_OBL = fitn(k, X_OBL, D);
            if Ffun_OBL < Ffun(i)
                X(i, :) = X_OBL;
                Ffun(i) = Ffun_OBL;
                if Ffun_OBL < Best_FF
                    Best_FF = Ffun_OBL;
                    Best_P = X_OBL;
                end
            end
     
             end

            
    %Update the convergence curve

    
    
        Conv_curve(C_Iter)=Best_FF;
    

    C_Iter=C_Iter+1;  % incremental iteration
   
end



