%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.39 null steering
% (a) Initial sinc-pattern
% (b) Pattern with a null of zero order imposed at u=0.22
% (c) With null of first order too
% (d) With null of second order too
% 
% Lillian Xiaolan Xu 3/24/99
% Last updated by K. Bell 7/22/01, 9/30/01, 10/17/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=27; 
us=0;     
u=0:0.001:1;
u1=0:0.001:1;
n=conj(-(N-1)/2:(N-1)/2)';
vs=exp(i*n.*pi.*us); 

% ----------------null constraints
nulls=[ 0.269 0.47 0.6971];   
vnull=exp(i*n.*pi.*nulls(1)); 
C=vnull;
if size(nulls,2)>1
   for m=2:size(nulls,2)
      C=[C exp(i*n.*pi.*nulls(m))]; 
   end;
end;

C1=[C i*n.*vnull];
C2=[C1 -n.*n.*vnull];

Pc=C*inv(C'*C)*C';
Pc1=C1*inv(C1'*C1)*C1';
Pc2=C2*inv(C2'*C2)*C2';

Wd=1/N*ones(N,1);
Wo=(Wd'*(eye(N)-Pc))';
Wo1=(Wd'*(eye(N)-Pc1))';
Wo2=(Wd'*(eye(N)-Pc2))';
b=0;     
for m=1:N
   b=b+Wd(m)*exp(-i*(-(N+1)/2+m)*pi*u1);
end;
bdb1=20*log10(abs(b));
b1=real(b);
b=0;     
for m=1:N
   b=b+Wo(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb2=20*log10(abs(b));
b2=real(b);
b=0;     
 
w3=Taylor(N,35,8);
w33=Taylor(N,36,8);
w3_1=(w33'*(eye(N)-Pc))';
for m=1:N
   b=b+w3(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
beam = abs(b)/max(abs(b));
bdb3=20*log10(beam);
b=0; 
for m=1:N
   b=b+w3_1(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
beam_1 = abs(b)/max(abs(b));
bdb3_1=20*log10(beam_1);
     

figure (1)
% subplot(2,2,1)
plot(u1,bdb1,'b:',u,bdb3,'r--',u,bdb3_1,'g-','MarkerSize',4,'LineWidth',2);
% plot(nulls(1),-80:10,'.')
hold on
plot(nulls(1),-80:10,'.','LineWidth',2);
plot(nulls(2),-80:10,'.','LineWidth',2);
plot(nulls(3),-80:10,'.','LineWidth',2);
axis([0 1 -80 10])
% title('(a) Uniform','Fontsize',14)
h1=xlabel('$u$','Fontsize',16);
h2=ylabel('Beam pattern (dB)','Fontsize',16);
h=legend('Uniform','Talor distribution','Nulls $u$=-0.43, 0.27, 0.70');
% line([nulls(1),nulls(1)],[-80 10]);
set(h1,'Interpreter','latex');
set(h,'Fontsize',14);
set(h,'Interpreter','latex');
grid on
hold on


% subplot(2,2,2)
% plot(u,bdb2,'LineWidth',2);
% hold on
% axis([-1 1 -80 10])


% plot(u,bdb3,'r-o','LineWidth',2);
% hold on
% axis([-1 1 -80 10])
% 
% plot(u,bdb3_1,'g->','LineWidth',2);
% hold on
% axis([-1 1 -80 10])
% 
% plot(nulls(1),-80:10,'.','LineWidth',2);
% plot(nulls(2),-80:10,'.','LineWidth',2);
% plot(nulls(3),-80:10,'.','LineWidth',2);
% title('(b) Zero-order null','Fontsize',14)
% h3=xlabel('$u$','Fontsize',16);
% h4=ylabel('Beam pattern (dB)','Fontsize',16);
% line([nulls(1),nulls(1)],[-80 10]);
% set(h3,'Interpreter','latex');
% h=legend('Uniform','Talor','Zero-order nulls $u$=-0.43 and 0.27');
% set(h,'Fontsize',14);
% set(h,'Interpreter','latex');
% grid on

% figure(3)
% % subplot(2,2,3)
% plot(u,bdb3);
% hold on
% axis([-1 1 -80 10])
% title('(c) First-order null','Fontsize',14)
% xlabel('\it u','Fontsize',14)
% ylabel('Beam pattern (dB)','Fontsize',14)
% line([nulls(1),nulls(1)],[-80 10]);
% grid on
% 
% figure(4)
% % subplot(2,2,4)
% plot(u,bdb4);
% hold on
% axis([-1 1 -80 10])
% title('(d) Second-order null','Fontsize',14)
% xlabel('\it u','Fontsize',14)
% ylabel('Beam pattern (dB)','Fontsize',14)
% line([nulls(1),nulls(1)],[-80 10]);
% grid on
