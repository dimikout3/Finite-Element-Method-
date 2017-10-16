function [ x ] = bicgstabPreco( A,b,max,M)
x=zeros(length(A),1);

r = b - A*x;
rhat=r;

rho=1;
alpha=1;
omega=1;

v=zeros(length(A),1);
p=zeros(length(A),1);

for iter= 1:max  
    rho2=rho;
    
    rho=rhat'*r;
    beta=( (rho/rho2)*(alpha/omega));
    p=r+beta*(p-omega*v);
    
    y=M\p;
    v=A*y;
    
    alpha=rho/(rhat'*v);
    
    s=r-alpha*v;
    
    z=M\s;
    t=A*z;
    
    w1=((M\t)'*(M\s));
    w2=(M\t)'*(M\t);
    omega=w1/w2;
    
    x=x+alpha*y+omega*z
    r=s-omega*t;
   
    if(norm(r)<1e-10*norm(rhat))   
        disp('Number of iterations is :');
         iter
        break;
    end    
end


