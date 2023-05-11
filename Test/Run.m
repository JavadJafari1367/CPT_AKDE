clc
close all
clear

data = xlsread('B&I CPT 2014.xlsx','Table 1');
data1=data([1:73 75:end],:);

Canliquefy=data1(:,3);
indexliq=find(Canliquefy==1);
indexnliq=find(Canliquefy==0);
indexmargin=find(Canliquefy==2);
qc1Ncs=data1([indexliq;indexnliq],22);
CSRM=data1([indexliq;indexnliq],27);

d = 2;

xi = qc1Ncs;
yi = CSRM;

MinX = min(xi);
MaxX = max(xi);
MinY = min(yi);
MaxY = max(yi);

gridX=MinX:1:MaxX;%
gridY=MinY:0.01:MaxY;
gridX = sort(gridX,'ascend');
gridY = sort(gridY,'ascend');

[B1 ind1 rem_ind1] = remove_similar_rows(gridX);
[B2 ind2 rem_ind2] = remove_similar_rows(gridY);

gridX = B1;
gridY = B2;

ci =[ones(1,180)';-ones(1,70)'];

XX =[xi yi];

figure(1)
plot(xi(ci==1),yi(ci==1),'.b')
hold on
plot(xi(ci==-1),yi(ci==-1),'.r')

X_1 =XX(ci==1,:);
N = size(X_1,1);
Sigma_Default = (max(X_1) - min(X_1)) .* ones(N,d);
Sigma_ML_H_d_1 = Sigma_Default;

X_0 =XX(ci==-1,:);
N = size(X_0,1);

Sigma_Default = (max(X_0) - min(X_0)) .* ones(N,d);
Sigma_ML_H_d_0 = Sigma_Default;

for it = 1:3
    MinX = min(xi);
    MaxX = max(xi);
    MinY = min(yi);
    MaxY = max(yi);
    
    Sigma_ML_H_d_1 =  BWE_Fast_ML(X_1,Sigma_ML_H_d_1,1);
    Sigma_ML_H_d_0 =  BWE_Fast_ML(X_0,Sigma_ML_H_d_0,1);
    
    KDE_1 = KDE2D(X_1, MinX , MaxX, MinY , MaxY, Sigma_ML_H_d_1);
    KDE_0 = KDE2D(X_0, MinX , MaxX, MinY , MaxY, Sigma_ML_H_d_0);
    
    [gridx22,gridx11] = meshgrid(gridY,gridX);
    
    figure
    [C,h]=contourf(gridx11,gridx22,KDE_1);
    clabel(C,h,'LabelSpacing',72,'FontSize',8)
    colormap(repmat((256:-1:0)'./256,1,3))
    
    hold on
    plot(xi(ci==1),yi(ci==1),'sr','MarkerSize',6)
    
    xlabel('q_c_1_N_c_s');
    ylabel('CSR_M_=_7_._5_,_σ_v_=_1_a_t_m')
    
    figure
    [C,h]=contourf(gridx11,gridx22,KDE_0);
    clabel(C,h,'LabelSpacing',72,'FontSize',8)
    colormap(repmat((256:-1:0)'./256,1,3))
    
    hold on
    plot(xi(ci==-1),yi(ci==-1),'ob','MarkerSize',6)
    
    xlabel('q_c_1_N_c_s');
    ylabel('CSR_M_=_7_._5_,_σ_v_=_1_a_t_m')
    
end
PL=KDE_1./(KDE_1+KDE_0);
figure
[C,h]=contourf(gridx11,gridx22,PL);
clabel(C,h,'LabelSpacing',72,'FontSize',8)
colormap(repmat((256:-1:0)'./256,1,3))
xlabel('q_c_1_N_c_s');
ylabel('CSR_M_=_7_._5_,_σ_v_=_1_a_t_m')
hold on
plot(xi(ci==1),yi(ci==1),'.r')
hold on
plot(xi(ci==-1),yi(ci==-1),'ob')