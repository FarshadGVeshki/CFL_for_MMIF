function [A1,A2] = SOMP_C(X1,X2,D1,D2,K,Eps)
% coupled SOMP (supp(A1) = supp(A2) with coupled dictionaries)
% Only atom slection is done jointly:  [~,ind] = max(abs(res1'*D1)+abs(res2'*D2));
% Coefficients are updated separately:   a1 = (D1(:,IND))\x1; and  a2 = (D2(:,IND))\x2;

X1 = single(X1);
X2 = single(X2);
D1 = single(D1);
D2 = single(D2);


[~, np] = size(X1);
[m, na] = size(D1);
EPS = Eps.^2;
A1 = single(zeros(na,np));
A2 = A1;
IND = zeros(K,1);
inds = [];

warning('off','all')

for i=1:np
    x1 = X1(:,i);
    x2 = X2(:,i);
    
    k = 0;
    res1 = x1;
    res2 = x2;
    res_norm1 = res1'*res1;
    res_norm2 = res2'*res2;
    a1 = [];
    a2 = [];
    while k<K && res_norm1>EPS && res_norm2>EPS
        k = k+1;
        [~,ind] = max(abs(res1'*D1)+abs(res2'*D2));
        IND(k) = ind;
        inds = IND(1:k);

        if k ==1
            a1 = D1(:,inds)'*x1;
            a2 = D2(:,inds)'*x2;
        else            
            a1 = D1(:,inds)\x1;
            a2 = D2(:,inds)\x2;
        end
        
        
        res1 = x1 - D1(:,inds)*a1;
        res2 = x2 - D2(:,inds)*a2;
        res_norm1 = res1'*res1;
        res_norm2 = res2'*res2;
    end
    
    if k>0
        A1(inds,i) = a1;
        A2(inds,i) = a2;
    end
    
end
warning('on','all')

A1 = double(A1);
A2 = double(A2);


