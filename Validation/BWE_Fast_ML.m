function  SigmaML = BWE_Fast_ML(X , Sigma, NN)
N = size(X,1);
D = size(X,2);

MinSigma = ((max(X) - min(X)) / N/10).^(1/2);

SigmaML = Sigma;
for it = 1 : NN
    M = zeros(N,1);
    for j = 1 : N
        M(j) =0 ;
        for i = 1 : N
            if(i ~= j)
                MS = 1;
                EXP = 0;
                for d2 = 1 : D
                    MS = MS * Sigma(i,d2);
                    EXP = EXP - ((X(i,d2)-X(j,d2)).^2 / (2 * Sigma(j,d2)^2));
                end
                
                M(j) = M(j) + ((1/(((pi .* 2)^(D/2)) .* MS)) .* exp(EXP));
            end
        end
        if(M(j) == 0)
            M(j) = 1;
        end
    end
    
    for k = 1 : N
        for d1 = 1 : D
            a = 0;
            c=0;
            
            for j = 1 : N
                MS = 1;
                EXP = 0;
                for d2 = 1 : D
                    MS = MS * Sigma(k,d2);
                    EXP = EXP - (((X(k,d2)-X(j,d2)).^2) / (2 * (Sigma(k,d2)^2)));
                end
                
                
                Wjk = ((1/(((pi .* 2)^(D/2)) .* MS)) .* exp(EXP))/M(j);
                
                if(k ~= j)
                    a = a +(Wjk * ((X(k,d1) - X(j,d1))^2));
                    c =  c + Wjk;
                end
                
            end
            
            if (a > 0) && (c ~= 0) && ~(isnan(a) || isnan(c))
                Sk =(a / c)^(0.5);
                if(Sk < MinSigma(d1))
                    Sk = MinSigma(d1);
                end
            else
                Sk = Sigma(k,d1);
            end
            
            SigmaML(k,d1) = Sk;
%             LEN = (Sigma(k,d1).*2);
%             
%             b = 0;
%             d = (LEN-(2*Sk))./(LEN - Sk);
%             a = -a;
%             r =roots([d c b a]);
%             r = r(imag(r)==0);
%             if i== j
%                 r = r(real(r)>0);
%             end
%             r = r(r< LEN);
%             
%             if(length(r) == 1)
%                 SigmaML(k,d1) =r(1);
%             elseif(length(r) > 1)
%                 SigmaML(k,d1) = mean(r);
%             else
%                 SigmaML(k,d1) = Sigma(k,d1);
%             end
%             
%             if(SigmaML(k,d1) < MinSigma(d1))
%                 SigmaML(k,d1) = MinSigma(d1);
%             end
        end
    end
    
    Sigma = SigmaML;
end
end
