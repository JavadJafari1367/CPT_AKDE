function  KDE = KDE2D(x , MinX , MaxX, MinY , MaxY, Sigma)
S = 140;

StepX = (MaxX-MinX)/S;
MinToMaxX = MinX : StepX : MaxX ;

StepY = (MaxY-MinY)/S;
MinToMaxY = MinY : StepY : MaxY ;

MinX = MinToMaxX(1);
MaxX = MinToMaxX(end);

MinY = MinToMaxY(1);
MaxY = MinToMaxY(end);

N = size(x,1);

KDE = zeros(length(MinToMaxX),length(MinToMaxY));

IndexX = 0;

P1 = (1./(((pi .* 2).^(3/2)) .* Sigma(:,1) .* Sigma(:,2) .* N));
P2 = (2 .*[Sigma(:,1).^2 Sigma(:,2).^2 ]);

for XX = MinToMaxX
    IndexX = IndexX + 1;
    IndexY = 0;
    for YY = MinToMaxY
        IndexY = IndexY + 1;
        
        K = [XX YY];
        
        KDE(IndexX,IndexY) =  sum(P1 .* exp(-1 .* sum((((x - K).^2) ./ P2),2)));
    end
end

KDE = KDE ./ sum(sum(KDE));

end



