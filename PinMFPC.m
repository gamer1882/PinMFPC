function tY= PinMFPC(X,iY,c1,c2,sigma,mt,tau)
% X is data,Each column is a sample,n*m
% Y is getted by Initialization.
% Y must be 1,2,...k.
% Cluster till convergence.
%totalw=zeros(n,num);
%mt the iter number of recursive
%num is the number of clustes.
n=size(X,1);%  
m=size(X,2);% 
flag=0;
alliter=0;
py=iY;%the lable of samples 
while flag==0 && alliter<100
    alliter=alliter+1;
    tY=py;
    L=unique(tY);
    num=length(L);
    wti=zeros(n,mt,num);
    CenterX=zeros(n,num);
    dis=zeros(m,num);%Initialize distance matrix from each sample to each cluster
    for i=1:num
        Ai=X(:,tY==L(i));%Ai is composed of samples belonging to the i-th cluster
        Bi=X(:,tY~=L(i));%Bi is composed of samples that do not belong to the i-th cluster
        CenterX(:,i)=mean(Ai,2); %average
        for t=1:mt
            if t==1
                wt=zeros(n,1);
            elseif norm(wt)~=0
                wt=wt/norm(wt);
            end
            Ai=Ai-(wt*wt')*Ai;
            Bi=Bi-(wt*wt')*Bi;
            w0=RPTW0(Ai); %Initialize w0
            wt=computewi(Ai,Bi,w0,c1,c2,sigma,tau);%wi
            wti(:,t,i)=wt;
        end   
    end
    for i=1:num
        M=wti(:,:,i)'*X-repmat(wti(:,:,i)'*CenterX(:,i),1,m);
        for j=1:m
            dis(j,i)=norm(M(:,j));%distance from each sample to each cluster
        end
    end 
    [~,py]=min(dis,[],2);%Update Y
    if getAC(tY,py)>0.9999
        flag=1;
    end
    %fprintf('.');
end
%fprintf('\n');
end
function ac=getAC(ty,py)
% ty: true class
% py: predict class
if size(ty,1)~=size(py,1) || size(ty,2)~=size(py,2)
    error('The size of two vector is not the same in getAC.m');
end
m=size(ty,1);
tM=zeros(m*(m-1)/2,1);
pM=zeros(m*(m-1)/2,1);
num=1;
column=1;
for i=1:m-1
    column=column+1;
    j=column;
    while j~=m+1 
        if ty(i,1)==ty(j,1)
            tM(num,1)=1;
        end
        if py(i,1)==py(j,1)
            pM(num,1)=1;
        end
        num=num+1;
        j=j+1;
    end
end
ac=size(find(tM==pM),1)/(m*(m-1)/2);
end

function w = RPTW0(A)%Initialize w0
meam_A = mean(A,2);
BarA = A - meam_A;
H = BarA * BarA';
[V,D] = eig(H);
[~,n] = min(diag(D));
w = real(V(:,n));
end