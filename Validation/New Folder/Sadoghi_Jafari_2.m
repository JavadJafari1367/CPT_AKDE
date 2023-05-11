clear
clc
load KDE

KDE = KDE_1 - KDE_0;
Sx = 200;
Sy = 100;
StepX = (MaxX-MinX)/Sx;
MinToMaxX = MinX - (StepX * 10) : StepX : MaxX + (StepX * 10);

StepY = (MaxY-MinY)/Sy;
MinToMaxY = MinY - (StepY * 10) : StepY : MaxY + (StepY * 10);

% StepZ = (MaxZ-MinZ)/Sz;
% MinToMaxZ = MinZ - (StepZ * 10) : StepZ : MaxZ + (StepZ * 10);

MinX = MinToMaxX(1);
MaxX = MinToMaxX(end);

MinY = MinToMaxY(1);
MaxY = MinToMaxY(end);

e = 0.000001;

R1 = zeros (1000000,2);
R1_Index = 1;

B1 = zeros (1000000,2);
B1_Index = 1;

R2 = zeros (1000000,2);
R2_Index = 1;

B2 = zeros (1000000,2);
B2_Index = 1;

for IndexX = 1 : size(KDE_1,1) - 1
    for IndexY = 1 : size(KDE_1,2) - 1
%         for IndexZ = 1 : size(KDE_1,3) - 1
            
            
            if (KDE_1(IndexX,IndexY) >= e)
                B1 (B1_Index,:) = [MinToMaxX(IndexX) MinToMaxY(IndexY)];
                B1_Index = B1_Index +1;
            end
            if (KDE_0(IndexX,IndexY) >= e)
                R1 (R1_Index,:) = [MinToMaxX(IndexX) MinToMaxY(IndexY)];
                R1_Index = R1_Index +1;
            end
            
            if (KDE_1(IndexX,IndexY) >= e) || (KDE_0(IndexX,IndexY) >= e)
                if (KDE(IndexX,IndexY) > 0)
                    B2 (B2_Index,:) = [MinToMaxX(IndexX) MinToMaxY(IndexY)];
                    B2_Index = B2_Index +1;
                end
                if (KDE(IndexX,IndexY) < 0)
                    R2 (R2_Index,:) = [MinToMaxX(IndexX) MinToMaxY(IndexY)];
                    R2_Index = R2_Index +1;
                end
                
            end
            
%         end
    end
end



figure(1)
clf
plot(xi(ci==1),yi(ci==1),'*b')
hold on
plot(xi(ci==-1),yi(ci==-1),'*r')
plot(R2(:,1),R2(:,2),'.','color',[1 0 0])
plot(B2(:,1),B2(:,2),'.','color',[0 0 1])




figure (2)
clf
plot(R1(:,1),R1(:,2),'.','color',[1 0 0])
figure (3)
clf
plot(B1(:,1),B1(:,2),'.','color',[0 0 1])

drawnow

