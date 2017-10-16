function [x]=CG(A,f,e,N)
x=zeros(size(f));
r0=f;
p=r0;

for i=1:N
    Ap=A*p;
    a=(r0'*r0)/(p'*Ap);
    x=x+a*p;
    r1=r0-a*Ap;
    if norm(r1)<e*norm(f)
        break;
    end
    b=(r1'*r1)/(r0'*r0);
    p=r1+b*p;
    r0=r1;
end

if i==N
    disp('Max Iterations Reached');
else
    disp(['Iterations: ' num2str(i) ' Residual: ' num2str(norm(r1))]);
end
end
