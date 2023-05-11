function  KDE = KDE2D(x , MinX , MaxX, MinY , MaxY, Sigma)

gridX=MinX:1:MaxX;
gridY=MinY:0.01:MaxY;
gridX = sort(gridX,'ascend');
gridY = sort(gridY,'ascend');

[B1 ind1 rem_ind1] = remove_similar_rows(gridX);
[B2 ind2 rem_ind2] = remove_similar_rows(gridY);

gridX = B1;
gridY = B2;

N = size(x,1);

KDE = zeros(length(gridX),length(gridY));

IndexX = 0;

P1 = (1./(((pi .* 2).^(2/2)) .* Sigma(:,1) .* Sigma(:,2)  .* N));
P2 = (2 .*[Sigma(:,1).^2 Sigma(:,2).^2]);

for XX = gridX
    IndexX = IndexX + 1;
    IndexY = 0;
    for YY = gridY
        IndexY = IndexY + 1; 
            K = [XX YY];          
            KDE(IndexX,IndexY) =  sum(P1 .* exp(-1 .* sum((((x - K).^2) ./ P2),2)));
    end
end

KDE = KDE ./ sum((sum(KDE)));

end



