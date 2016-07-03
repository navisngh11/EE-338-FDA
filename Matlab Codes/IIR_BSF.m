% Filter number 73
% Navjot Singh
% 130110071
% Chebyshev

% Global Parameters
delta=0.15;                
f_sample = 100000;
f_message = 45000;
t_width = 2000;

% Given analog frequencies corresponding to filter number
B_l= 16300;
B_h= 26300;
B_l_p= B_l - t_width;
B_h_p= B_h + t_width;
analog_specs= [B_l_p, B_l, B_h, B_h_p];  

% Conversion to digital specification
digital_specs= (analog_specs*2*pi)/f_sample;

% Bilinear transformation to analog domain
% In the form (Omega_p1,Omega_s1,Omega_s2,Omega_p2)
analog_freq= tan(digital_specs/2);   

% Calculation of paramteres
% For Bandstop to Lowpass transfromation
Omega_0= sqrt(analog_freq(1)*analog_freq(4));
B= analog_freq(4)-analog_freq(1);

% Bandpass to Lowpass transformation
% Using Omega_lp = B/((Omega_0^2/Omega_bp)- Omega_bp)
% Gives Passband at +-1
analog_lp_freq= B./(((Omega_0^2)./analog_freq)-analog_freq);
% Choosing the stringent of the two stopbands
Omega_lp_s= min(abs(analog_lp_freq(2)),abs(analog_lp_freq(3)));
Omega_lp_p= 1;

% Design of Chebyshev Filter
D1= 1/(1-delta)^2 -1;
D2= 1/delta^2 -1;
epsilon=sqrt(D1);
N_temp = (acosh(sqrt(D2)/epsilon))/(acosh(Omega_lp_s/Omega_lp_p));
N = ceil(N_temp);           %Order of Chebyshev filter
 
% Poles of Chebyshev filter (all poles)
% Are of the form s= j*Omega_p*cos(A_k+jB_k)
% Hence real part = Omega_p*sin(A_k)*sinh(B_k)
% Hence imaginary part = Omega_p*cos(A_k)*cosh(B_k)
% Where A_k=(2k+1/2N)*pi
% and B_k= asinh(1/epsilon)/N
poles_total_real=zeros(1,2*N);
poles_total_imag=zeros(1,2*N);
for k= 0:1:(2*N)-1
    A_k = (((2*k)+1)*pi)/(2*N);
    B_k = asinh(1/epsilon)/N;
    poles_total_real(k+1) = Omega_lp_p*sin(A_k)*sinh(B_k);
    poles_total_imag(k+1) = Omega_lp_p*cos(A_k)*cosh(B_k);
end
poles_total_real;
poles_total_imag;
poles_total=poles_total_real+j*poles_total_imag;

% Taking poles on the LHP ;
poles_lp=[];
for k= 1:(2*N)
     if(poles_total_real(k)<0)
         poles_lp=[poles_lp,poles_total_real(k)+j*poles_total_imag(k)];
     end
end
poles_lp;
% [z,p,k]=cheb1ap(4,1.4114); To check poles (Ripple = 10log(1+epsilon^2))


% Preparing Transfer Function
% Numerator depends upon parity on N
% To get the gain at Omega_p equal to 1
syms s 
if(mod(N,2)==0)
    num_gain=1/sqrt(1+epsilon^2);
else
    num_gain=1;
end
num_lp(s)= (-1^N)*num_gain*prod(poles_lp)*(s^0);
denom_lp(s)=s^0;
real(poly(poles_lp));     %return polynomial coeffs whose roots are poles_lp
for k=1:N                 %for checking
   denom_lp(s)= (s-poles_lp(k))*denom_lp(s);
end
sym2poly(num_lp(s));
sym2poly(denom_lp(s));
% Taking only real part
real(sym2poly(denom_lp(s)));
num_lp(s);
denom_lp(s);
T_analog_lp(s)=num_lp(s)/denom_lp(s);


% % Lowpass to Bandstop Transformation
syms s_bs;
s= (B*s_bs)/(s_bs^2 + Omega_0^2);
T_analog_bs(s_bs) = T_analog_lp(s);
% Finding the numerator and denominator polynomials
[num_bs,denom_bs]=numden(T_analog_bs);
real(sym2poly(num_bs));       
real(sym2poly(denom_bs));


% Converting to Digital domain
% Let the variable be z^-1
% Hence s_bs -> (1-z^-1/1+z^-1)
% P.S : coeffs start from z^-n
syms invz
s_bs = (1-invz)/(1+invz);
T_digital_bs(invz)=T_analog_bs(s_bs);
[num,denom]=numden(T_digital_bs);
coeff_num=real(sym2poly(num)) ;  
coeff_den=real(sym2poly(denom));
freqz(coeff_num,coeff_den);
digital_specs/pi

% This part is to convert the scale to linear from dB
% Courtesy : MathWorks FAQs
ax=get(gcf,'Children'); 			
li=get(ax(1),'Children'); 
r=get(ax(1),'YLabel'); 
set(r,'String','Magnitude') 
ydata=get(li,'Ydata'); 
y=10.^(ydata/20); 
set(li,'Ydata',y);