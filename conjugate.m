function bestx =conjugate(H,c2,A,K,tau)
% min 0.5x'Hx+Kx+c2e' max(e-Ax,-tau(e-Ax))
tol=1e-3;
n=size(H,1);
x0=ones(n,1);
g0=grad(H,c2,A,K',x0,tau);
d0=-g0;
if norm(g0)<tol
    bestx=x0;
%    bestf=obj(H,c2,A,K,x0,tau);
    return;
end
x0_1=x0;
dH=d0'*H;
dHd=dH*d0;
b=(dH*x0+K*d0)/dHd;
c=2*c2/dHd;
q=A*x0;
s=A*d0;
a0=Smart(b,c,q,s,tau);
g0_1=g0;
x0=x0+a0*d0;
g0=grad(H,c2,A,K',x0,tau);
ite=1;
while norm(x0-x0_1)>tol && ite<500
    x0_1=x0;
    beta=g0'*(g0-g0_1)/(g0_1'*g0_1);
    d0=-g0+beta*d0;
    dH=d0'*H;
    dHd=dH*d0;
    if dHd<1e-6
        break;
    end
    b=(dH*x0+K*d0)/dHd;
    c=2*c2/dHd;
    q=A*x0;
    s=A*d0;
    ind=find(s~=0);
    a0=Smart(b,c,q(ind,1),s(ind,1),tau);
    g0_1=g0;
    x0=x0+a0*d0;
    g0=grad(H,c2,A,K',x0,tau);
    ite=ite+1;
end
bestx=x0;
%bestf=obj(H,c,A,K,bestx,tau);
end

function val=obj(H,c2,A,K,x0,tau)
val=0.5*x0'*H*x0+c2*sum(max((1-A*x0),-tau*(1-A*x0)))+K*x0;
end

function right=grad(H,c2,A,K,x,tau)
tol=1e-6;
% grad(x)=Hx+K'-c2A'sign(e-Ax)
right=zeros(size(A,1),1);
Ax=A*x;
Axe=1-Ax;
right(Axe<-tol,1)=-tau;
right(Axe>tol,1)=1;
right=H*x-c2*A'*right+K;
end