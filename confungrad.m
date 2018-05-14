function [c ceq DC DCeq] = confungrad(x)

global H Sigma b a x0bar N numcon


k=numcon;
m=2*N;

eps = x(1:k);
u = x(k+1:k+m);
c=zeros(k,1);
z=zeros(k,1);
Sig1=zeros(k,1);
for i = 1:k
    z(i) = b(i)-a(:,i)'*(x0bar + H*u);
    Sig1(i) = a(:,i)'*Sigma*a(:,i);
    c(i) = 1 - eps(i) - normcdf(z(i)/sqrt(Sig1(i)));
end
ceq = [];

if nargout >2
    DC = zeros(k+m,k);
    DC(1:k,:) = -1.*eye(k);
     for i = 1:k
        d1(i) = exp(-z(i)^2/(2*Sig1(i)));
        d2(i) = sqrt(2*pi*Sig1(i));
        DC(k+1:k+m,i) = (d1(i)/d2(i))*H'*a(:,i);
     end  
    
    DCeq=[];
end
        