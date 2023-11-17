% ------------------------TAYLOR SERIES -----------------------
% -------------------------------------------------------------------------
%written by Thirusabaresaan.P,VIT UNIVERSITY.
% mail to: thirusabari2000@gmail.com
% -------------------------------------------------------------------------
clear all
clc
syms x y z
f= input('enter function  ') %enter the functions
x0=input('initial points  ') %enter the limits
m=input('enter number of variables  ') %number of variables 1,2,3
%single variable
if m==1
    n=input('enter how many terms  ') %enter required number of terms for taylor series
    f_0=subs(f,x,x0)
for i=1:n 
    if f~=0
    f1(i)=diff(f,x)
    f11(i)=subs(f1(i),x,x0)
    B(i)=(x-x0)^(i)
    f=f1(i);
    end
end
a1=transpose(f11);
f_1=B*a1
f_ty=f_0+f_1 %taylor series non simplifed
f_taylor=(expand(f_ty)) %taylor series simplifed
% two variables    
elseif m==2
gf=gradient(f) %gradient
sub_gf=subs(gf,[x,y],[x0(1) x0(2)])
a=[x-x0(1) y-x0(2)] %value of 'X'
f0=subs(f,[x,y],[x0(1) x0(2)])
f1=a*sub_gf
hf=hessian(f) %hessian
sub_hf=subs(hf,[x,y],[x0(1) x0(2)])
f21=a*sub_hf
b=transpose(a)
f2=(1/2)*f21*b
f_taylor=simplify(f0+f1+f2) %taylor series
% three variables
elseif m==3
gf=gradient(f)%gradient
sub_gf=subs(gf,[x,y,z],[x0(1) x0(2) x0(3)])
a=[x-x0(1) y-x0(2) z-x0(3)] %value of 'X'
f0=subs(f,[x,y,z],[x0(1) x0(2) x0(3)])
f1=a*sub_gf
hf=hessian(f) %hessian
sub_hf=subs(hf,[x,y,z],[x0(1) x0(2) x0(3)])
f21=a*sub_hf
b=transpose(a)
f2=(1/2)*f21*b
f_Taylor=simplify(f0+f1+f2) %taylor series
end