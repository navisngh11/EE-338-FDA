% Filter number 73
% Navjot Singh
% 130110071
% Butterworth

% Global Parameters
delta=0.15;                
f_sample = 100000;
f_message = 45000;
t_width = 2000;

% Given analog frequencies corresponding to filter number
B_l= 14900;
B_h= 24900;
B_l_s= B_l - t_width;
B_h_s= B_h + t_width;
analog_specs= [B_l_s, B_l, B_h, B_h_s]; 

% Conversion to digital specification
digital_specs= (analog_specs*2*pi)/f_sample;

% Bilinear transformation to analog domain
% In the form (Omega_s1,Omega_p1,Omega_p2,Omega_s2)
analog_freq= tan(digital_specs/2);   

% Calculation of paramteres
% For Bandpass to Lowpass transfromation
Omega_0= sqrt(analog_freq(2)*analog_freq(3));
B= analog_freq(3)-analog_freq(2);

% Bandpass to Lowpass transformation
% Using Omega_lp = (Omega_bp/B) - (Omega_0^2)/(B*Omega_bp)
% Gives Passband at +-1
analog_lp_freq= (analog_freq/B) - ((Omega_0^2)./(B*analog_freq));
% Choosing the stringent of the two stopbands
Omega_lp_s= min(abs(analog_lp_freq(1)),abs(analog_lp_freq(4)));
Omega_lp_p= 1;

% Design on Butterworth Filter
D1= 1/(1-delta)^2 -1;
D2= 1/delta^2 -1;
N_temp = 0.5*log(D2/D1)/log(Omega_lp_s/Omega_lp_p);
N = ceil(N_temp);          %Order of Butterworth filter
omega_t1=Omega_lp_p/(D1^(1/(2*N)));
omega_t2=Omega_lp_s/(D2^(1/(2*N)));
Omega_c=(omega_t1+omega_t2)/2;

% Poles of Butterworth filter (all poles)
% Are of the form s=j*Omega_c*exp(j(2k+1)pi/2N)
% Hence real part = -Omega_c*sin((2k+1)pi/2N)
% Hence imaginary part = Omega_c*cos((2k+1)pi/2N)
poles_total_real=zeros(1,2*N);
poles_total_imag=zeros(1,2*N);
for k= 0:(2*N)-1
    poles_total_real(k+1) = -1*Omega_c*sin((((2*k)+1)*pi)/(2*N));
    poles_total_imag(k+1) = Omega_c*cos((((2*k)+1)*pi)/(2*N));
end

% Taking poles on the LHP Omega_c^N;
poles_lp=[];
for k= 1:(2*N)
     if(poles_total_real(k)<0)
         poles_lp=[poles_lp,poles_total_real(k)+j*poles_total_imag(k)];
     end
end


% Preparing Transfer Function
syms s 
num_lp(s)=Omega_c^N * (s^0);
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


% % Lowpass to Bandpass Transformation
syms s_bp;
s= (s_bp^2 + Omega_0^2)/(B*s_bp);
T_analog_bp(s_bp) = T_analog_lp(s);
% Finding the numerator and denominator polynomials
[num_bp,denom_bp]=numden(T_analog_bp);
real(sym2poly(num_bp));       
real(sym2poly(denom_bp));


% Converting to Digital domain
% Let the variable be z^-1
% Hence s_bp -> (1-z^-1/1+z^-1)
% P.S : coeffs start from z^-n
syms invz
s_bp = (1-invz)/(1+invz);
T_digital_bp(invz)=T_analog_bp(s_bp);
[num,denom]=numden(T_digital_bp);
coeff_num=real(sym2poly(num));   
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