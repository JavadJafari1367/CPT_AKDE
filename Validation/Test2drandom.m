MinX = -7;
MaxX = 7;
MinY = -7;
MaxY = 7;


x1 = 0
y1 = 0
s11 = 1
s12 = 1

x2 = 4
y2 = 2
s21 = 0.2
s22 = 1

x3 = -2
y3 = -4
s31 = 1
s32 = 0.2

R1 = mvnrnd([x1 y1],[s11 0; 0 s12],100);
R2 = mvnrnd([x2 y2],[s21 0; 0 s22],100);
R3 = mvnrnd([x3 y3],[s31 0; 0 s32],100);

X = [R1 ;R2 ;R3];

d = 2;
N = length(X);

Sigma_Default = ones(N,d) .* [max(X) - min(X)];
Sigma_ML_H_Best = Sigma_Default;

Sigma_ML_H_Best =  BWE_Fast_ML(X,Sigma_ML_H_Best,8);
KDE_1 = KDE2D(X, MinX , MaxX, MinY , MaxY,Sigma_ML_H_Best);
KDE_2 = KDE2D(X, MinX , MaxX, MinY , MaxY,ones(N,d) .* [0.1 1]);
KDE_3 = KDE2D(X, MinX , MaxX, MinY , MaxY,ones(N,d) .* [1 0.1]);
KDE_4 = KDE2D(X, MinX , MaxX, MinY , MaxY,ones(N,d) .* [0.7 0.7]);



MinToMaxX = -7:0.1:7;

MinToMaxY = -7:0.1:7;

KDE = zeros(length(MinToMaxX),length(MinToMaxY));

IndexX = 0;
for XX = MinToMaxX
    IndexX = IndexX + 1;
    IndexY = 0;
    for YY = MinToMaxY
        IndexY = IndexY + 1;
        
        K1 = [XX YY] - [x1 y1];
        K2 = [XX YY] - [x2 y2];
        K3 = [XX YY] - [x3 y3];
        
        KDE(IndexX,IndexY) =                      (1/(((pi .* 2).^(2/2)) .* s11 .* s12)) .* exp(-1 .* sum(((K1.^2) ./ (2 .*[s11^2 s12^2]))));
        KDE(IndexX,IndexY) = KDE(IndexX,IndexY) + (1/(((pi .* 2).^(2/2)) .* s21 .* s22)) .* exp(-1 .* sum(((K2.^2) ./ (2 .*[s21^2 s22^2]))));
        KDE(IndexX,IndexY) = KDE(IndexX,IndexY) + (1/(((pi .* 2).^(2/2)) .* s31 .* s32)) .* exp(-1 .* sum(((K3.^2) ./ (2 .*[s31^2 s32^2]))));
    end
end

sum(sum(KDE))

KDE = KDE ./ sum(sum(KDE));


figure (1)
clf
hold on
contour(MinToMaxX,MinToMaxY,KDE',20,'r')
contour(MinToMaxX,MinToMaxY,KDE_1',20,'b')
legend("True Density","Adaptive Bandwidth")
xlabel("X")
ylabel("Y")

figure (2)
clf
hold on
contour(MinToMaxX,MinToMaxY,KDE',20,'r')
contour(MinToMaxX,MinToMaxY,KDE_2',20,'b')
legend("True Density","Bandwidth = (1.0,0.1)")
xlabel("X")
ylabel("Y")
figure (3)
clf
hold on
contour(MinToMaxX,MinToMaxY,KDE',20,'r')
contour(MinToMaxX,MinToMaxY,KDE_3',20,'b')
legend("True Density","Bandwidth = (0.1,1.0)")
xlabel("X")
ylabel("Y")
figure (4)
clf
hold on
contour(MinToMaxX,MinToMaxY,KDE',20,'r')
contour(MinToMaxX,MinToMaxY,KDE_4',20,'b')
legend("True Density","Bandwidth = (1.0,1.0)")
drawnow
xlabel("X")
ylabel("Y")

% figure (5)
% hold on
% mesh(MinToMaxX,MinToMaxY,KDE'-KDE_1')
