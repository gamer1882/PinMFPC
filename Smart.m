function bestx=Smart(b,c,q,s,tau)
% One dimension searching:
% min x^2+2bx+ce' max(e-(sx+q),-tau e-(sx+q))
% where b,c are numbers and d, e, f are vectors and e is ones.
% **Note* f must have no 0, otherwise, set d(f~=0,1) and f(f~=0,1) as input 
% **Note* length(f) should be lower than 50000, otherwise may have problem.
val=(1-q)./s; 
k=0; % coef of tmpx in
n=length(s);
% maxval=max(max(abs(s)),max(abs(q)));
% maxval=max(max(abs(b),abs(c)),maxval);
% tol=maxval*1e-8;
tol=1e-3;
for i=1:n
    if s(i)>0 && -b>val(i) || s(i)<0 && -b<val(i)
        k=k+tau*s(i);
    elseif s(i)>0 && -b<val(i) || s(i)<0 && -b>val(i)
        k=k-s(i);
    end
end
if abs(k)<tol
    bestx=-b;
    return;
else
    coefl=zeros(n,1);
    coefr=zeros(n,1);
    for i=1:n
        if s(i)>0
            coefl(i)=-s(i);
            coefr(i)=tau*s(i);
        elseif s(i)<0
            coefl(i)=tau*s(i);
            coefr(i)=-s(i);
        end
    end
end
region_coef=zeros(n,1);
last_point=[];
if k<0
    [val,ind]=sort(val);
    coefl=coefl(ind,1);
    coefr=coefr(ind,1);
    region_coef(1)=sum(coefl);
    for i=2:length(region_coef) 
        new_coef=region_coef(i-1,1)-coefl(i-1)+coefr(i-1);
        if new_coef<0
            region_coef(i,1)=new_coef;
        else
            last_point=i-1; %last point in pos for negative coef
            break;
        end
    end
    if isempty(last_point)
        last_point=n;
    end
    [~,right]=min(abs(abs(region_coef(1:last_point,1)-k)));%right=find(abs(region_coef(1:last_point,1)-k)<tol,1);
    mys=-(b+c/2*region_coef(right,1));
    while 1        
        if mys<=val(right,1) 
            bestx=mys;
            return;
        elseif right==last_point
            bestx=val(right,1);
            return;            
        else
            right=right+1;
        end
        mys=-(b+c/2*region_coef(right,1));
        if mys<=val(right-1,1)
            bestx=val(right-1,1);
            return;
        end
    end           
else
    [val,ind]=sort(val,'descend');
    coefl=coefl(ind,1);
    coefr=coefr(ind,1);
    region_coef(1)=sum(coefr);
    for i=2:length(region_coef) %number of region
        new_coef=region_coef(i-1,1)-coefr(i-1)+coefl(i-1);
        if new_coef>0
            region_coef(i,1)=new_coef;
        else
            last_point=i-1; %last point in pos for negative coef
            break;
        end
    end
    if isempty(last_point)
        last_point=n;
    end
    [~,left]=min(abs(abs(region_coef(1:last_point,1)-k)));%left=find(abs(region_coef(1:last_point,1)-k)<tol,1);
    mys=-(b+c/2*region_coef(left,1));
    while 1        
        if mys>=val(left,1)
            bestx=mys;
            return;
        elseif left==last_point
            bestx=val(left,1);
            return;            
        else
            left=left+1;
        end
        mys=-(b+c/2*region_coef(left,1));
        if mys>=val(left-1,1)
            bestx=val(left-1,1);
            return;
        end        
    end  
end

end