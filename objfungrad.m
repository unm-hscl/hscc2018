function [f dF] = objfungrad(x)
global N numcon

eps = x(1:numcon);
u = x(numcon+1:length(x));
f = sum(eps);

if nargout > 1
    dF=zeros(length(x),1);
    for i = 1:length(eps)
        dF(i)=1;
    end
    for i= length(eps)+1:length(u)+length(eps)
        dF(i)=0;
    end
end