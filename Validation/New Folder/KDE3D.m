function  KDE = KDE3D(x , MinX , MaxX, MinY , MaxY, MinZ , MaxZ, Sigma)
Sx = 50;
Sy = 100;
Sz = 50;
StepX = (MaxX-MinX)/Sx;
MinToMaxX = MinX - (StepX * 10) : StepX : MaxX + (StepX * 10);

StepY = (MaxY-MinY)/Sy;
MinToMaxY = MinY - (StepY * 10) : StepY : MaxY + (StepY * 10);

StepZ = (MaxZ-MinZ)/Sz;
MinToMaxZ = MinZ - (StepZ * 10) : StepZ : MaxZ + (StepZ * 10);

MinX = MinToMaxX(1);
MaxX = MinToMaxX(end);

MinY = MinToMaxY(1);
MaxY = MinToMaxY(end);

MinZ = MinToMaxZ(1);
MaxZ = MinToMaxZ(end);

N = size(x,1);

KDE = zeros(length(MinToMaxX),length(MinToMaxY),length(MinToMaxZ));

IndexX = 0;

P1 = (1./(((pi .* 2).^(3/2)) .* Sigma(:,1) .* Sigma(:,2)  .* Sigma(:,3) .* N));
P2 = (2 .*[Sigma(:,1).^2 Sigma(:,2).^2 Sigma(:,3).^2]);

for XX = MinToMaxX
    IndexX = IndexX + 1;
    IndexY = 0;
    for YY = MinToMaxY
        IndexY = IndexY + 1;
        IndexZ = 0;
        for ZZ = MinToMaxZ
            IndexZ = IndexZ + 1;
            
            K = [XX YY ZZ];
            
            KDE(IndexX,IndexY,IndexZ) =  sum(P1 .* exp(-1 .* sum((((x - K).^2) ./ P2),2)));
        end
    end
end

KDE = KDE ./ sum(sum(sum(KDE)));

end



