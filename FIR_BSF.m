% Filter number 73
% Navjot Singh
% 130110071
% FIR Filter using Kaiser Window

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

% Parameters for Kaiser Window
A = -1*20*log10(delta)
omega_t= min(digital_specs(2)-digital_specs(1),digital_specs(4)-digital_specs(3));
N_1 = (A-8)/(2*2.285*omega_t);
N = ceil(N_1)+6;   % Adjusted to get overshoot within Tolerance
if A>50
    alpha=0.1102*(A-8.7);
else
    if(A>=21)
        alpha=0.5842*((A-21)^0.4)+(0.07886*(A-21));
    else
        alpha=0;
    end
end
beta= alpha/N;

% Preparing Ideal Impulse Response
omega_c1=0.5*(digital_specs(1)+digital_specs(2));
omega_c2=0.5*(digital_specs(3)+digital_specs(4));
h_ideal= zeros((2*N)+1,1);
for k=-N:1:N
    h_ideal(k+N+1)=-1*(sin(omega_c2*k)-sin(omega_c1*k))/(pi*k);
end
h_ideal(N+1)=1-(omega_c2-omega_c1)/pi;
h_ideal;
% Generate Kaiser window coefficients
k_window=kaiser((2*N)+1,alpha);

% Windowing the Ideal Impulse Response
h_filter=h_ideal.*k_window;
freqz(h_filter);
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