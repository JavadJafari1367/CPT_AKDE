function [B ind rem_ind] = remove_similar_rows(A)
temp=[]; c = 2;
ind=[];rem_ind=[];
temp(1,:) = A(1,:);
ind(1) = 1;
for i=2:size(A,1)
    flag = 0;
    for j=1:size(temp,1)
        if(A(i,:) == temp(j,:))
            flag = 1;
            rem_ind = [rem_ind; [i j]];
            break;
        end
    end
    if(flag==0)
        temp(c,:) = A(i,:);
        ind(c) = i;
        c = c+1;
    end
end
B = temp;