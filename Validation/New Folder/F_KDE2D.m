function  KDE = F_KDE2D(x,Sigma ,K)
N = size(x,1);
P1 = (1./(((pi .* 2).^(2/2)) .* Sigma(:,1) .* Sigma(:,2) .* N));
P2 = (2 .*[Sigma(:,1).^2 Sigma(:,2).^2]);
KDE = sum(P1 .* exp(-1 .* sum((((x - K).^2) ./ P2),2)));
end
