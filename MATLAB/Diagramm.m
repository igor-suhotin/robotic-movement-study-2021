classdef Diagramm
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Parameter;
        Points;
        ClassOfPoints;
        kLength;
        muLength;
%         BoundsPoints;
%         BoundsCounts;
    end
    
    methods
        function D = Diagramm(parameter)
            
%             BoundsPoints = cell(4);
%             BoundsPoints{1}(1,1)=1;
%             BoundsPoints{1}(1,2)=2;
%             BoundsPoints{1}(2,1)=1;
%             BoundsPoints{1}(2,2)=2;
%             finish;
            D.Parameter=parameter;
            D.Points=[];
            D.ClassOfPoints=[];
            D.kLength=0;
            D.muLength=0;
            kNum=51; % 51;
            muNum=38; % 38;
%             for k=0:0.02:1
%                 for mu=1:0.08:4  %mu=parameter:0.08:3+parameter
            for I=1:kNum
                k=(I-1)*(1.0/(kNum-1));
                for J=1:muNum  %mu=parameter:0.08:3+parameter
                    mu=1+(J-1)*(3.0/(muNum-1));
                    fprintf('length(D.Points)= %d ', length(D.Points));
                    D.Points((I-1)*muNum+J,1)=k;
                    fprintf(' %d ', length(D.Points));
                    D.Points((I-1)*muNum+J,2)=mu;
                    fprintf(' %d \n', length(D.Points));
                    D.ClassOfPoints((I-1)*muNum+J,1)=-1;
                    D.ClassOfPoints((I-1)*muNum+J,2)=-1;
%                     if D.muLength==0
%                         muNum=muNum+1;
%                     end
                end
%                 kNum=kNum+1;
                if D.muLength==0
                    D.muLength=muNum;
                end
            end
            D.kLength=kNum;
            fprintf('length(D.Points)= %d \n', length(D.Points));
%             finish;
        end
        
        function D = OnlyOneClassify(D,I)
            S=-1;
            N=-1;
            k=D.Points(I,1);
            mu=D.Points(I,2);
            [S, N] = count_of_stops_and_napravlenie(k, mu, D.Parameter);
            
%             if k<0.2
%                 S=0;
%             elseif k<0.5
%                 S=1;
%             elseif k<0.8
%                 S=2;
%             elseif k<1
%                 S=3;
%             end
%             if mu<1+D.Parameter
%                 N=0;
%             elseif mu<2+D.Parameter
%                 N=1;
%             elseif mu<3+D.Parameter
%                 N=2;
%             end
            
            D.ClassOfPoints(I,1)=S;
            D.ClassOfPoints(I,2)=N;
        end
        
        function D = AllClassify(D)
            for I=1:length(D.Points)
                D=OnlyOneClassify(D,I);
            end
        end

        function D = krest(k, myu)
            plot(k, myu, '+', 'Color', [0 0 0], 'LineWidth', 1);
            fprintf('krest  k=%f  myu=%f \n', k, myu);
        end
        

        function Postroit(D)
            disp('Postroit');
            index=@(kNum,muNum) (kNum-1)*D.muLength+muNum;
            figure
            xlim([0 1])
            ylim([1 4])
            hold on
            title(['b=',num2str(D.Parameter)]);
            xlabel('k') 
            ylabel('mu')

            FaceColors = [ [1 0 0]; [0 1 0]; [0 0 1]; [1, 165.0/255, 0]; [1 1 1]]; %red, green, blue, orange, white
            ConturColors = [ [1 0 1]; [1 1 0]; [0 1 1] ]; %magenta (розовый), yellow, cyan (сине-зелёный)

%             D.BoundsPoints = [ ];            
%             BoundsCounts = [ ];
%             BoundsPoints = cell(length(FaceColors));
%             for I=1:length(FaceColors)
%                 BoundsPoints(I) = {[]};
%                 BoundsCounts(I) = 0;
%             end

            [bp1,bp2,bp3,bp4] = deal([],[],[],[]);
            [bc1,bc2,bc3,bc4] = deal(0);

%             for I=1:length(D.Points)
            for kNum=1:D.kLength
              for muNum=1:D.muLength
                I = index(kNum,muNum);
                FaceColor = FaceColors(D.ClassOfPoints(I,1)+1,:);
                ConturColor = ConturColors(D.ClassOfPoints(I,2)+2,:);
%                if D.ClassOfPoints(I,1) == 0
%                    FaceColor=[1 0 0]; %red
%                elseif D.ClassOfPoints(I,1) == 1
%                    FaceColor=[0 1 0]; %green
%                elseif D.ClassOfPoints(I,1) == 2
%                    FaceColor=[0 0 1]; %blue
%                elseif D.ClassOfPoints(I,1) == 3
%                    FaceColor=[1 1 1]; %white
%                elseif D.ClassOfPoints(I,1) == -1
%                    FaceColor=[0 0 0]; %black
%                end
%                if D.ClassOfPoints(I,2) == -1
%                    ConturColor=[1 0 1];  %magenta розовый
%                elseif D.ClassOfPoints(I,2) == 0
%                    ConturColor=[1 1 0]; %yellow
%                elseif D.ClassOfPoints(I,2) == 1
%                    ConturColor=[0 1 1];  %cyan сине-зелёный
%                elseif D.ClassOfPoints(I,2) == 2
%                    ConturColor=[0 0 0];  %black
%                end
                plot(D.Points(I,1),D.Points(I,2),'o','Color',ConturColor,'MarkerFaceColor',FaceColor,'LineWidth',2);
%                 fprintf('I=%d\n', I);
                leftIndex=index(kNum-1,muNum);
                if abs(D.Points(I,1)-0.850000) < 0.0001 && abs(D.Points(I,2)-1.52941) < 0.0001
                    disp(D.Points(I,2));
                end
                
                
%                 if kNum>1 && D.ClassOfPoints(I,1) == D.ClassOfPoints(leftIndex,1) || ...
%                    kNum<D.kLength && D.ClassOfPoints(I,1) == D.ClassOfPoints(index(kNum+1,muNum),1) || ...
%                    muNum>1 && D.ClassOfPoints(I,1) == D.ClassOfPoints(I-1,1) || ...
%                    muNum<D.muLength && D.ClassOfPoints(I,1) == D.ClassOfPoints(I+1,1)
                if not ((kNum<D.kLength && D.ClassOfPoints(index(kNum+1,muNum),1) == 4) || (muNum<D.muLength && D.ClassOfPoints(I+1,1) == 4))
                    dist=1;
                    if D.ClassOfPoints(I,1) == 3
                        dist = 2;
                    end
                    if kNum>1 && D.ClassOfPoints(I,1)~=4 && D.ClassOfPoints(I,1) > D.ClassOfPoints(leftIndex,1) && D.ClassOfPoints(I,1) <= dist+D.ClassOfPoints(leftIndex,1)
                        P = PointDihotomy(D.Parameter, D.Points(leftIndex,1), D.Points(leftIndex,2), D.ClassOfPoints(leftIndex,:), D.Points(I,1), D.Points(I,2), D.ClassOfPoints(I,:), 15);
%                         plot(P(1), P(2), '+', 'Color', [1 0 0], 'LineWidth', 1);
                        if D.ClassOfPoints(I,1) == 0
                            bc1 = bc1+1;
                            [bp1(bc1, 1), bp1(bc1, 2)] = deal(P(1), P(2));
                        elseif D.ClassOfPoints(I,1) == 1
                            bc2 = bc2+1;
                            [bp2(bc2, 1), bp2(bc2, 2)] = deal(P(1), P(2));
                        elseif D.ClassOfPoints(I,1) == 2
                            bc3 = bc3+1;
                            [bp3(bc3, 1), bp3(bc3, 2)] = deal(P(1), P(2));
                        elseif D.ClassOfPoints(I,1) == 3
                            bc4 = bc4+1;
                            [bp4(bc4, 1), bp4(bc4, 2)] = deal(P(1), P(2));
                        end
    %                    fprintf('  krest  k=%f  myu=%f P(2)=%f\n', P(1), D.Points(I,2), P(2));
    %                     krest(P(1), D.Points(I,2));
                        if false && muNum>1 && D.ClassOfPoints(leftIndex,1) == D.ClassOfPoints(leftIndex-1,1)
                            P = PointDihotomy(D.Parameter, D.Points(leftIndex-1,1), D.Points(leftIndex-1,2), D.ClassOfPoints(leftIndex-1,:), D.Points(I,1), D.Points(I,2), D.ClassOfPoints(I,:), 16);
%                             plot(P(1), P(2), '+', 'Color', [1 0 0], 'LineWidth', 1);
                            if D.ClassOfPoints(I,1) == 0
                                bc1 = bc1+1;
                                [bp1(bc1, 1), bp1(bc1, 2)] = deal(P(1), P(2));
                            elseif D.ClassOfPoints(I,1) == 1
                                bc2 = bc2+1;
                                [bp2(bc2, 1), bp2(bc2, 2)] = deal(P(1), P(2));
                            elseif D.ClassOfPoints(I,1) == 2
                                bc3 = bc3+1;
                                [bp3(bc3, 1), bp3(bc3, 2)] = deal(P(1), P(2));
                            elseif D.ClassOfPoints(I,1) == 3
                                bc4 = bc4+1;
                                [bp4(bc4, 1), bp4(bc4, 2)] = deal(P(1), P(2));
                            end
    %                         BoundsCounts(D.ClassOfPoints(I,1)+1)=BoundsCounts(D.ClassOfPoints(I,1)+1)+1;
    %                         BoundsPoints{D.ClassOfPoints(I,1)+1}(BoundsCounts(D.ClassOfPoints(I,1)+1))=P;
    %                         D.BoundsPoints(D.BoundsCounts(D.ClassOfPoints(I,1)+1),2)=D.Points(I,2);
                        end
                    end
                

                    % PointDihotomy( b, k_begin, myu_begin, C_begin, k_end, myu_end, C_end, steps)
                    if muNum>1 && D.ClassOfPoints(I,1)~=4 && D.ClassOfPoints(I,1) > D.ClassOfPoints(I-1,1) && D.ClassOfPoints(I,1) <= dist+D.ClassOfPoints(I-1,1)
                        P = PointDihotomy(D.Parameter, D.Points(I-1,1), D.Points(I-1,2), D.ClassOfPoints(I-1,:), D.Points(I,1), D.Points(I,2), D.ClassOfPoints(I,:), 15);
%                         plot(P(1), P(2), '+', 'Color', [0 0 0], 'LineWidth', 1);
                        if D.ClassOfPoints(I,1) == 0
                            bc1 = bc1+1;
                            [bp1(bc1, 1), bp1(bc1, 2)] = deal(P(1), P(2));
                        elseif D.ClassOfPoints(I,1) == 1
                            bc2 = bc2+1;
                            [bp2(bc2, 1), bp2(bc2, 2)] = deal(P(1), P(2));
                        elseif D.ClassOfPoints(I,1) == 2
                            bc3 = bc3+1;
                            [bp3(bc3, 1), bp3(bc3, 2)] = deal(P(1), P(2));
                        elseif D.ClassOfPoints(I,1) == 3
                            bc4 = bc4+1;
                            [bp4(bc4, 1), bp4(bc4, 2)] = deal(P(1), P(2));
                        end
%                    fprintf('krest  k=%f  myu=%f  P(1)=%f\n', D.Points(I,1), P(2), P(1));
%                    krest(D.Points(I,1), P(2));
                    end
                end
              end
            end
            
            if bc2 > 1
                bp2 = sortrows(bp2,1);
                plot2 = plot(bp2(:,1), bp2(:,2),'-','LineWidth',2);

            end
            if bc3 > 1 
                bp3 = sortrows(bp3,1);
                plot3 = plot(bp3(:,1),bp3(:,2),'-','LineWidth',2);
            end
            if bc4 > 1 
                bp4 = sortrows(bp4,1);
                plot4 = plot(bp4(:,1),bp4(:,2),'-','LineWidth',2);
                
            end
            
            disp('====== bp1 ======');
            disp(bp1);
            disp('====== bp2 ======');
            disp(bp2);
            disp('====== bp3 ======');
            disp(bp3);
            disp('====== bp4 ======');
            disp(bp4);
            
            k_on_mu = @(mu) 1/(D.Parameter*sqrt(mu^2-1));
            mu_on_k = @(k) sqrt(1+1/(D.Parameter*k)^2);
            
            x=k_on_mu(4):0.001:1;
            y=1:length(x);
            for I=1:length(x)
                y(I)=mu_on_k(x(I));
            end
            plot(x,y,'-','Color', [0 0 0],'LineWidth',2);
        end
        
    end
end

