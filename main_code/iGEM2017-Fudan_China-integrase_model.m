pb=0.01;
int=1;
pbi=0.0001;
lri1=0.0001;
lri2=0.0001;
lr=0.0001;
intA=0;
t=0:0.01:10.05;
dt=0.01;
beta=1;
alpha=log(2)/0.5;

for i=1:1:1005
%     if i>200 & beta>0.05
%         beta=beta-0.05
% 
%     end
    
    y=kinetics(pb(i),pbi(i),lri1(i),lri2(i),lr(i),int(i),intA(i),dt,beta,alpha)
    pb(i)=y(1);
    pbi(i)=y(2);
    lri1(i)=y(3);
    lri2(i)=y(4);
    lr(i)=y(5);
    int(i)=y(6);
    intA(i)=y(7);
    
    y=thermonew(pb(i),pbi(i),lri1(i),lri2(i),lr(i),int(i))
    pb(i+1)=y(1);
    pbi(i+1)=y(2);
    lri1(i+1)=y(3);
    lri2(i+1)=y(4);
    lr(i+1)=y(5);
    int(i+1)=y(6);
    intA(i+1)=intA(i);
    
%     y=thermo(pb(i),pbi(i),lri1(i),lri2(i),lr(i),int(i))
%     pb(i+1)=y(1)
%     pbi(i+1)=y(2)
%     lri1(i+1)=y(3)
%     lri2(i+1)=y(4)
%     lr(i+1)=y(5)
%     int(i+1)=y(6)
%     intA(i+1)=intA(i)

end
figure
hold on

plot(t,(0.0101-lri2)/0.0001,'LineWidth',2)
%plot(t,pbi)
plot(t,(lri2-0.0001)/0.0001,'LineWidth',2)
%plot(t,(lri2+lri1+lr)/0.0001,'LineWidth',2)
legend('pb','lr')
xlabel('time,h')
ylabel('proportion,%')
axis([0,10,0,101]);
%plot(t,lr)
%plot(t,int)
%plot(t,lr+lri1+lri2)
%plot(t,intA)
grid on

function y=kinetics(pb,pbi,lri1,lri2,lr,int,intA,dt,beta,alpha)
kpr=6
kmr=2.14
kpsyn=0.006
kmsyn=0.017
kbi=0.0001
lrtot=lr+lri1+lri2
pbtot=pb+pbi
dnatot=lrtot+pbtot
a=pbi-kpr*pbi*dt+kmr*lri1*dt-alpha*dt*pbi
b=lri1-(kpsyn+kmr)*lri1*dt+kmsyn*lri2*dt+kpr*pbi*dt-alpha*dt*lri1
c=lri2+kpsyn*lri1*dt-kmsyn*lri2*dt-alpha*dt*lri2
d=int+beta*dt-alpha*int*dt
e=intA+beta*dt*pbtot/dnatot-alpha*intA*dt
f=pb+pbi*alpha*dt
h=lr+lri1*alpha*dt+lri2*alpha*dt
pbi=a
lri1=b
lri2=c
int=d
intA=e
pb=f
lr=h

y=[pb,pbi,lri1,lri2,lr,int,intA]

end

function y=thermonew(pb,pbi,lri1,lri2,lr,int)
kbi=0.0001
z(1)=(pb*int^4-kbi*pbi)/(kbi+int^4)
z(2)=(lr*int^4-kbi*lri2)/(kbi+int^4)
pb=pb-z(1)
pbi=pbi+z(1)
lr=lr-z(2)
lri2=lri2+z(2)
y=[pb,pbi,lri1,lri2,lr,int]
end



function y=thermo(pb,pbi,lri1,lri2,lr,int)
kbi=0.0001
A=[1/pb+1/pbi+16/int,16/int;16/int,1/lr+1/lri2+16/int]
b=[log(kbi*pbi)-log(pb*int^4);log(kbi*lri2)-log(lr*int^4)]
z=A\b

pbs=pb+z(1)
pbis=pbi-z(1)

lri2s=lri2-z(2)
lrs=lr+z(2)

ints=int+4*z(1)+4*z(2)
if pbs>0 & pbis>0 & lri2s>0 & lrs>0 & ints>0
    pb=pbs
    pbi=pbis
    lri2=lri2s
    lr=lrs
    int=ints
else
    A=1/pb+1/pbi+16/int
    b=log(pbi*int^4)-log(kbi*pb)
    z=b/A
    pbs=pb+z
    pbis=pbi-z
    ints= int-4*z
if pbs>0 & pbis>0 & ints>0
    pb=pbs
    pbi=pbis
    int=ints

else
    A=1/lr+1/lri2+16/int
    b=log(kbi*lri2)-log(lr*int^4)
    z=b/A
    
    lri2s=lri2-z
    lrs=lr+z
    ints=int+4*z
    if lri2s>0 & lrs>0 &ints>0
        lri2=lri2s
    lr=lrs
    int=ints
    end
    
end
end


y=[pb,pbi,lri1,lri2,lr,int]
end

