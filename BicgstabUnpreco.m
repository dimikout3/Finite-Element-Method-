function [ x ] = BicgstabUnpreco( A,b,max)
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
    
    v=A*p;
    alpha=rho/(rhat'*v);
    
    s=r-alpha*v;

    t=A*s;
   
    omega=(t'*s)/(t'*t);
    
    x=x+alpha*p+omega*s;
    r=s-omega*t;
   
    if(norm(r)<1e-10*norm(rhat))   
        disp('Number of iterations is :');
         iter
        break;
    end    
end


