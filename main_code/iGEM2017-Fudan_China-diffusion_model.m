r=0:0.01:1;
dr=r(2)-r(1);
t=0:0.0001:10;
dt=t(2)-t(1);
c=zeros(1,101);
alpha=log(2)/0.5;
beta=0.1;
d=0.0278;  %΢��ƽ��/s
for j=1:100000
   
    c=fickslaw(c,dr,alpha,beta,d,r,dt);
end
c=c/c(1)
plot(r,c,'LineWidth',2)
xlabel('relative position')
ylabel('relative concentration')
legend('c')
axis([0,1.01,0,1.01])
text(0.17,c(18),'\leftarrow average distance(0.17,0.6577)')
grid on
    function c=fickslaw(cc,dr,alpha,beta,d,r,dt)
    dcdr=zeros(1,101);
    d2cdr=zeros(1,100);
    ctemp=zeros(1,101);
    cc(1)=cc(1)+beta*dt;%/(4*pi*r(2)*r(2)*r(2)/3);       %��Դ��λʱ��������µķ��ӣ�ϵ��beta
    detacc=0;
    sumctemp=0;
    sumcc=0;
    
for i=2:101
   
dcdr(i)=(cc(i-1)-cc(i))/dr;
end

for i=2:100
    
d2cdr(i)=(dcdr(i)-dcdr(i+1))/dr;

end

for i=2:100
ctemp(i)=cc(i)+(d*d2cdr(i)+2*d*dcdr(i)/r(i))*dt;   %��ɢ����
%sumctemp=sumctemp+ctemp(i)*r(i)*r(i)*dr*4*pi ;  %�µ�Ũ�����
%sumcc=sumcc+cc(i)*r(i)*r(i)*dr*4*pi;            %ԭŨ�����
        

end
detacc=sumctemp-sumcc;                  %����һ��λʱ���е�Դ�������ɢ��Ũ��
ctemp(1)=cc(1);          %-ctemp(2)/; %��Դ��������ɢ���µ�Ũ�ȱ仯
ctemp=ctemp-alpha*ctemp*dt;              %ͳһ���⣬ϵ��alpha
c=ctemp;                                %����Դ�����Ũ�ȷ��ص�c��
    end
    
    
    
