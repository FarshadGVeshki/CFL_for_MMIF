function D = DCT(K,bb)    

    Pn=ceil(sqrt(K));
    D=zeros(bb,Pn);
    for k=0:1:Pn-1
        V=cos([0:1:bb-1]'*k*pi/Pn);
        if k>0, V=V-mean(V); end    
    
        D(:,k+1)=V/norm(V);
    end
    D=kron(D,D);
    
    D = D(:,1:K);