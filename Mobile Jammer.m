%%%%%%%%%%%%% The main.m file %%%%%%%%%%%%%%%
clear;
% Parameters of the Gaussian filter:
n1=10;sigma1=3;n2=10;sigma2=3;theta1=0;
% The amplitude of the noise:
noise=0.1;
[w,map]=gifread('lena.gif');
x=ind2gray(w,map);
filter1=d2gauss(n1,sigma1,n2,sigma2,theta);
x_rand=noise*randn(size(x));
y=x+x_rand;
f1=conv2(x,filter1,'same');
rf1=conv2(y,filter1,'same');
figure(1);
subplot(2,2,1);imagesc(x);
subplot(2,2,2);imagesc(y);
subplot(2,2,3);imagesc(f1);
subplot(2,2,4);imagesc(rf1);
colormap(gray);
%%%%%%%%%%%%%% End of the main.m file %%%%%%%%%%%%%%%

%Initialization D=1 (D = interleaving depth)
clear;
Ts=1/1e3; % tx rate = 10k bps
NA=10000;
A=floor(rand(1,NA)+0.5);
%A=randint(1,NA);
input=A;
G=[1 0 0; 1 0 1; 1 1 1];
k=1;

N=length(channel_output);
SNRdB= -10:5:20;
SNR=10.^(SNRdB/10);
Counter=zeros(1,length(SNRdB));
loop=0;

%************************************************* **************
% AWGN + fading *
%************************************************* **************

Fd=66; % Actual fading rate
fm=round(Fd*NA*Ts);

%Generate two inputs of 1 x N matrix
xin=randn(1,(NA+2)*3);
yin=randn(1,(NA+2)*3);

%Loop for all values of frequencies and truncate at just before f=fm
for f=1fm-1)
filter(f)=1.5/(pi*fm*sqrt(1-(f/fm)^2));
end
for f=(fm-1)NA+2)*3
filter(f)=0;
end

%Inputs of the filter is multiplied with the filter transfer function
xout=xin.*filter;
yout=yin.*filter;

%Perform inverse Fourier Transform on both output
xifft=real(ifft(xout));
yifft=real(ifft(yout));

% normaliation so that E{R^2}= 1 / in Watt = 0 dB
std_x=std(xifft);
std_y=std(yifft);
xifft=xifft/std_x/sqrt(2);
yifft=yifft/std_y/sqrt(2);

%Sum the square the of the IFFT and take the square root of it
rx=sqrt(xifft.^2+yifft.^2);

%************************************************* ***************
%Calculate for all value of SNRdB in AWGN + fading
%************************************************* ***************
loop=0;
for b=sqrt(SNR) %b is the signal level
loop=loop+1;
dfad=0;
for n=1NA+2)*3 %n = no of bits
x=randn;
t=channel_output(n);
if t>=0.5
a=b*rx(n);
r=a+x;
if r<0
dfad=dfad+1; % error has occured, increase counter
Df(n)=0;
else
Df(n)=1;
end
else
a=(-b)*rx(n);
r=a+x;
if r>0
dfad=dfad+1;
Df(n)=1;
else
Df(n)=0;
end
end
Counterf(loop)=dfad; % w/o decode
end
output=viterbi(G,k,Df); % E is output of decoder
zf=0;
v=1;
for v=v:1:length(A);
if xor(input(1,v),output(1,v))==1
zf=zf+1;
end
end
Counter_wcf(loop)=zf;
end% loop for signal (b) in AWGN + fading

%plot simulation in AWGN + fading without coding
Pef=Counterf/N;
%semilogy(SNRdB,Pef,'b');
%hold on;
%plot simulation in AWGN + fading with coding
Pe1f=Counter_wcf/NA;
save D1_66SNRdBPe1fPef
%semilogy(SNRdB,Pe1f,'b'); %D=1 no interleaving
%hold on;

clear;
Ts=1/1e3; % tx rate = 10k bps
NA=10000;
A=floor(rand(1,NA)+0.5);
%A=randint(1,NA);
input=A;
G=[1 0 0; 1 0 1; 1 1 1];
k=1;

N=length(channel_output);
SNRdB= -10:5:20;
SNR=10.^(SNRdB/10);
Counter=zeros(1,length(SNRdB));
loop=0;