function [ x ] = cgunprecon( A,b,iter,M )
x=zeros(length(A),1);
r=b-A*x;
r0=r;
z=M\r;
p=z;
for i=1:iter
    rpast=r;
    zpast=z;
    a1=(r'*z);
    a2=(p'*A*p);
    a=a1/a2;
    x=x+a*p;
    r=r-a*A*p;
    if norm(r)<(1e-10*norm(r0)) 
        disp('number of iterations is :')
        i
        break
    end
    z=M\r;
    b=(z'*r)/(zpast'*rpast);
    p=z+b*p;
end
if i==iter 
    disp('Max iterations reached at :')
    iter
end

