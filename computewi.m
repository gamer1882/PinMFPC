function wi = computewi(Ai,Bi,w0,c1,c2,sigma,tau)
% Ai is the sample data belonging to the i-th cluster
% c1,c2,sigma and tau are paramater
% w0 is the initial value
X=[Ai,Bi];%all sample data
[n,mAi]=size(Ai);
[~,mBi]=size(Bi);
eAi=ones(mAi,1);
iter=0;
tol=0.001;
B=(Bi-repmat(mean(Ai,2),1,mBi))';
HA=Ai-mean(Ai,2)*eAi';
H=eye(n)+HA*HA'*c1;  
XX=X*X';
while iter<100
    iter=iter+1;
    Gi=diag(sign((B*w0)));
    Hi=sign(w0'*XX*w0-1);%in artical is hi(Lowercase h).   
    K=sigma*Hi*w0'*XX;
    A=Gi*B;
    wi =conjugate(H,c2,A,K,tau);
    som=norm(wi-w0);%wi-w0
    if som<tol
        break;
    end
    w0=wi;
end
%iter
end