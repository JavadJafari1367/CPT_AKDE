% close all
clear;
clc;

% dataT = xlsread('New Moss data.xlsx','Sheet1');
% data = dataT(2:end,1:3);
% datatest = dataT(2:183,1:3);
% rng(1,'twister');
Zhoa=xlsread('Zhao2022database','NCEER1998');
Canliquefy=Zhoa(:,1);%
indexliq=find(Canliquefy==1);
indexnliq=find(Canliquefy==0);
M=Zhoa([indexliq;indexnliq],2); % moment magnitude (M)
amax=Zhoa([indexliq;indexnliq],3);% peak ground acceleration (g)
depth=Zhoa([indexliq;indexnliq],4);% Critical depth range (m)
gwt=Zhoa([indexliq;indexnliq],5); % ground water table (m)
sigmavo=Zhoa([indexliq;indexnliq],6); % vertical total stress (kPa)
sigmapvo=Zhoa([indexliq;indexnliq],7); % vertical effective stress (kPa)
qc=Zhoa([indexliq;indexnliq],8); % cone tip resistance (kPa)
fs=Zhoa([indexliq;indexnliq],9); % sleeve friction (kPa)
pa=101.325;% 1 atm kPa
n=ones(1,length(fs)).';
for i=1:length(fs)
    deltan=1;
    while abs(deltan)>=0.01
        Cn=(pa/sigmapvo(i,1))^n(i,1);
        Q=((qc(i,1)-sigmavo(i,1))/pa)*Cn;
        F=(100*fs(i,1))/(qc(i,1)-sigmavo(i,1));
        Ic=sqrt((3.47-log10(Q))^2+(1.22+log10(F))^2);
        %         nj=0.50.*(Ic<=1.64)+((Ic-1.64)*0.30+0.50).*(Ic>1.64 & Ic<3.30)+1*(Ic>=3.30);
        nj=(0.381*Ic+0.05*(sigmapvo(i,1)/pa)-0.15);nj=nj*(nj<=1)+1*(nj>1);
        deltan=nj-n(i,1);
        n(i,1)=nj;
    end
end
Cn=(pa./sigmapvo).^n;Cn=Cn.*(Cn<=1.7)+1.7.*(Cn>1.7);
Q=((qc-sigmavo)/pa).*Cn;
F=(100*fs)./(qc-sigmavo);
Ic=sqrt(((3.47-log10(Q)).^2+(1.22+log10(F)).^2));
for i=1:length(fs)
    if Ic(i,1)<2.50
        Kc(i,1)=1*(Ic(i,1)<=1.64)+(-0.403*Ic(i,1).^4+5.581*Ic(i,1).^3-21.63*Ic(i,1).^2+33.75*Ic(i,1)-17.88).*(Ic(i,1)>1.64 & Ic(i,1)<=2.50)+1*((Ic(i,1)>1.64 & Ic(i,1)<2.36) & F(i,1)<0.5);
    end
    if Ic(i,1)>2.50 && Ic(i,1)<2.70
        Kc(i,1)=(6*10^-7)*(Ic(i,1)).^16.76;
    end
end
qc1Ncs=Kc.*Q;
% CRR=(93*(qc1Ncs/1000).^3+0.08).*(qc1Ncs>=50 & qc1Ncs<=160)+(0.833*(qc1Ncs/1000)+0.05).*(qc1Ncs<50);
rqc1Ncs=qc1Ncs.*(qc1Ncs<=211)+211.*(qc1Ncs>211);
% rqc1Ncs=qc1Ncs.*(qc1Ncs<=200)+211.*(qc1Ncs>200);
Dr=-85+76*log10(rqc1Ncs);
Csigma=1./(37.3-8.27*(rqc1Ncs).^0.264);Csigma=Csigma.*(Csigma<=0.3)+0.30.*(Csigma>0.30);
Ksigma=1-Csigma.*log(sigmapvo/pa);Ksigma=Ksigma.*(Ksigma<=1.1)+1.1.*(Ksigma>1.1);
MSFmax=1.09+(qc1Ncs/180).^3; MSFmax=MSFmax.*(MSFmax<=2.2)+2.2.*(MSFmax>2.2);
MSF=1+(MSFmax-1).*(8.64*exp(-M/4)-1.325);
az=-1.012-1.126*sin((depth/11.73)+5.133);
bz=0.106+0.118*sin((depth/11.28)+5.142);
rd=exp(az+bz.*M);
CSR=0.65*(sigmavo./sigmapvo).*amax.*rd;
CSRM=CSR./(MSF.*Ksigma);

d = 2;

xi = Q;%data(:,1);% qc1
yi = CSR;%data(:,2);% CSR


MinX = min(xi);
MaxX = max(xi);
MinY = min(yi);
MaxY = max(yi);


% ci =data(:,3);
ci =[ones(1,194)';-ones(1,101)'];


XX_All =[xi yi];
clf
hold on

IDX = randperm(size(XX_All,1));

Train = fix(size(XX_All,1) * .80);
XX_Train = XX_All(IDX(1:Train),:);
XX_Test = XX_All(IDX(Train +1 :end),:);

C_Train = ci(IDX(1:Train),:);
C_Test = ci(IDX(Train +1 :end),:);

% figure(1)
% plot3(XX_Train(C_Train==1,1),XX_Train(C_Train==1,2),XX_Train(C_Train==1,3),'.b')
% hold on
% plot3(XX_Train(C_Train==0,1),XX_Train(C_Train==0,2),XX_Train(C_Train==0,3),'.r')

X_1 =XX_Train(C_Train==1,:);
N = size(X_1,1);
Sigma_Default = (max(X_1) - min(X_1)) .* ones(N,d);
Sigma_ML_H_d_1 = Sigma_Default;

X_0 =XX_Train(C_Train==-1,:);
N = size(X_0,1);
Sigma_Default = (max(X_0) - min(X_0)) .* ones(N,d);
Sigma_ML_H_d_0 = Sigma_Default;
deltap=[];Aucp=deltap;Accp=deltap;MCCp=deltap;
threshold=[];jpoop=[];
for it = 2
    Sigma_ML_H_d_1 =  BWE_Fast_ML(X_1,Sigma_ML_H_d_1,it);
    Sigma_ML_H_d_0 =  BWE_Fast_ML(X_0,Sigma_ML_H_d_0,it);
    jp=1;
    dj=0.005;
    for j = 1:dj:100
        thr = j/100;
        
        TP = 0;
        TN = 0;
        FP = 0;
        FN = 0;
        
        for i = 1 : size(XX_Test,1)
            F_C0 = F_KDE2D(X_0,Sigma_ML_H_d_0,XX_Test(i,:));
            F_C1 = F_KDE2D(X_1,Sigma_ML_H_d_1,XX_Test(i,:));
            P1 = F_C1 / (F_C0+F_C1);
            %                         P1 = (194/295)*F_C1 / ((101/295)*F_C0 + (194/295)*F_C1);
            if(P1>=thr) && C_Test(i) == 1
                TP = TP + 1;
            elseif(P1<thr) && C_Test(i) == -1
                TN = TN + 1;
            elseif(P1>=thr) && C_Test(i) == -1
                FP = FP + 1;
            elseif(P1<thr) && C_Test(i) == 1
                FN = FN + 1;
            end
        end
        Sen(it,jp) = TP / (TP + FN);
        Spe(it,jp) = TN / (TN + FP);
        Acc(it,jp) = (TP + TN) / (TP +FP + TN +FN);
        Pre_P(it,jp) = TP / (TP + FP);
        Pre_N(it,jp) = TN / (TN + FN);
        MCC(it,jp) = (TP*TN-FN*FP)/sqrt((TP+FN)*(TN+FP)*(TP+FP)*(TN+FN));
        
        F_Me_P(it,jp) = (2 * Pre_P(it,jp) * Sen(it,jp)) / (Pre_P(it,jp)+ Sen(it,jp));
        F_Me_N(it,jp) = (2 * Pre_N(it,jp) * Sen(it,jp)) / (Pre_N(it,jp)+ Sen(it,jp));
        
        
        Rec(it,jp) = TP/(TP + FN);
        FPR(it,jp) = FP/(FP + TN);
        delta=FPR(it,jp)-1+Rec(it,jp);
        if abs(delta)<0.075
            jpoop=[jpoop;jp];
            deltap=[deltap;delta];
            threshold=[threshold;j];
%             disp([delta j])
        end
%         plot(FPR(it,jp),Rec(it,jp),'.b')
        %     plot(j,Sen,'.b')
        %     plot(j,Spe,'.k')
        %     plot(j,Pre_P,'.r')
        %     plot(j,Pre_N,'.y')
        %     plot(j,MCC,'.g')
        
        jp=jp+1;
    end
    disp(it)
    Aucp(it)=-trapz(FPR(it,:),Rec(it,:));
    a=find(abs(deltap)==min(abs(deltap)));
    if ~(isempty(a))
    indoop=jpoop(a);%+a(1);%+round(a(end)-a(1))/2;
%     indoop=round(indoop(1));
    oop=roundn((indoop(1)-1)*dj+1,-2);
%     indoop=round(((roundn(threshold(index),-2)-1)/dj)+1);
    disp([Acc(it,indoop(1)) MCC(it,indoop(1)) Aucp(it) oop])
    Accp=[Accp;Acc(it,indoop(1))];
    MCCp=[MCCp;MCC(it,indoop(1))];
    end
    
    %     r1(it) = (TP)/(TP+FN);
    %     r2(it) = (FP)/(TP+FN);
    %     r3(it) = (TN)/(TN+FP);
    %     r4(it) = (FN)/(TN+FP);
    %
    %     d1(it) = (TP);
    %     d2(it) = (FP);
    %     d3(it) = (TN);
    %     d4(it) = (FN);
end

% figure (1)
% clf
% hold on
% plot(r1)
% plot(r2)
% plot(r3)
% plot(r4)
% legend("TP","FP","TN","FN")
%
% figure (2)
% clf
% hold on
% plot(d1)
% plot(d2)
% plot(d3)
% plot(d4)
% legend("TP","FP","TN","FN")

x=FPR(2,:);
y=Rec(2,:);
z=1:dj:100;zp=z/100;
surf([x(:) x(:)], [y(:) y(:)], [zp(:) zp(:)], ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', 2);            % Make a thicker line
view(2);   % Default 2-D view
colorbar;  % Add a colorbar
colormap(jet)