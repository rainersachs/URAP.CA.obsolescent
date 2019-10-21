r1=.5; r2 =1-r1; alph=3; A=.15
d=.001*1:1500
f=exp(r1*d)+A*exp(alph*r2*d)-1-A
g=exp(d)-1
h=A*(exp(alph*d)-1)
plot(d,f,type='l')
lines(d,g,col='red')
lines(d,h,col='green')
