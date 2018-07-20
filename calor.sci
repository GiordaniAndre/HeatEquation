//U_t = U_xx*eta(x) + U_x*eta'(x) + f(x,t)

function z=c_inicial(x)

z = x*(1-x)

endfunction

function z=funcao(x,t)
    
    z = x*(1-x)*exp(-t)
        
endfunction




function z=eta1(x)

    z = (x-1/2)^2

endfunction

function z=eta_x(x)

    z = 2*(x-1/2)

endfunction

function y=heat(L,tf,N_x,N_t,theta)

h=L/(N_x-1)
dt=tf/(N_t-1)


x=[0:h:L]   //meshgrid

a=ones(N_x,1)
b=ones(N_x,1)
c=ones(N_x,1)

    
    for k=1:N_x
    
        a(k) = ((-dt/(h^2))*(1-theta)*eta1(x(k)) + eta_x(x(k))*(dt/(2*h))*(1-theta))
        b(k) = (2*(dt/(h^2))*(1-theta)*eta1(x(k)))+1
        c(k) = ((-dt/(h^2))*(1-theta)*eta1(x(k)) - eta_x(x(k))*(dt/(2*h))*(1-theta))
    
    
    end

    a(1)   = 0 //doesn't matter the value, it's out of the matrix
    a(N_x) = 0 //0 for dirichlet
    
    c(1)   = 0 //0 for dirichlet
    c(N_x) = 0 //doesn't matter the value, it's out of the matrix
    
    b(1)   = 1
    b(N_x) = 1

    U0 = zeros(N_x,1)
    d = zeros(N_x,1)


    d(1)=1
    d(N_x)=-1
    
    
    for k=1:N_x
        U0(k)=c_inicial(x(k))        
    end
 
      for j=1:N_t
            for k=2:N_x-1
                d(k) = U0(k)
                d(k) = d(k) + (theta)*(U0(k+1) - 2*U0(k) + U0(k-1))*(dt/h^2)*eta1(x(k)) //segunda derivada
                d(k) = d(k) + (theta)*(U0(k+1) - U0(k-1))*eta_x(x(k))*(dt/(2*h)) //primeira derivada
                d(k) = d(k) + (theta)*dt*funcao(x(k),(j-1)*dt) //função f(x,t)
                d(k) = d(k) + (1-theta)*dt*funcao(x(k), (j)*dt) // função em t+1
            end
          
        
        U0     = TDMA(a,b,c,d)
        
          //  plot(x',U0)


        if modulo(j,int(N_t/100))==1
            plot(x',U0)
        end
    end

   
    y      = [x' U0]
    y      = [[0:1:N_x-1]' a b c d U0]
         plot(x',U0)
    y=[a b c d]
endfunction

function x=TDMA(a,b,c,d)

n=size(a,1) // recover the system order 

cl=zeros(n,1) //Inicialize cl  
dl=zeros(n,1) //Inicialize dl  
x=zeros(n,1)  //Inicialize x  

cl(1)=c(1)/b(1)  
for i=2:n-1  
    cl(i)=c(i)/(b(i)-a(i)*cl(i-1))  
end  

dl(1)=d(1)/b(1)  
for i=2:n  
    dl(i)=(d(i)-a(i)*dl(i-1))/(b(i)-a(i)*cl(i-1))  
end  

x(n)=dl(n)  
for i=n-1:-1:1  
    x(i)=dl(i)-cl(i)*x(i+1)  
end  

endfunction 

function tempo(ti,tf)
    for k=ti:.05:tf
        heat(1,k,101,401,.5)
    end
endfunction
