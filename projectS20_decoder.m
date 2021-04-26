clc, clear 
SNR = [0,2,4,6,8,10];
L = 14700;
N= 10;
%%%%%%%%%% C Matrix Generation %%%%%%%%%
c= [1, 0.3,0.1,0.01];
c= transpose(c);
C = convmtx(c,N);
C(N+1,:) =[];
C(N+2,:) = [];
C(N+1,:) =[];

%%%%%%%%%%% Zero filtering %%%%%%%
Wn = inv(C);
Wn_ZF = Wn(:,1); 
for a= 1:length(SNR)
%%%%%%%% A Generations %%%%%%%%%%
A= sqrt((10^(SNR(a)/10))/2);

%%%% Transmiting Signal X matrix Generation %%%%%
x_randvalu = randi(4,[N*L,1]);
QPSK_Set = [complex(A,A), complex(A,-A), complex(-A,A),complex(-A,-A)];
x= QPSK_Set(x_randvalu);
x = transpose(x);

%%%%% Complex Gaussian Noise V matrix generation %%%%
v = 1/sqrt(2)*(randn(N*L,1)+1i*randn(N*L,1));


%%%%%% ISI Channel Generation %%%%%
t_matrix = conv(x,c);
Y = t_matrix(1:N*L)+v;

%%%%%%%%%%%%%%%%%%%%%% ZERO- FILTERING %%%%%%%%%%%%%%%%%%%%%%%%

%%%% Applying ZF for Y generation for ZERO- FILTERING %%%%
Y_dash = conv(Y,Wn_ZF); 

%%%% Finding Min A, Closest X to Y'(n) for ZERO- FILTERING %%%%%%%%%
QPSK_Set = [complex(A,A), complex(A,-A), complex(-A,A),complex(-A,-A)];
for k = 1:N*L
    for q = 1:4
    
       
        all_Diff = Y_dash(k) - QPSK_Set(q);
        norm_diff = norm(all_Diff,2);
        norm_Store(q) = norm_diff;
    

    end
   [minvalue , minindx] = min(norm_Store);
    value(k) = minvalue;
    indx(k)= minindx; 
    
end
%%%%%%% X_Hat generation for ZERO- FILTERING %%%%%%%%%%%
X_hat = QPSK_Set(indx);
X_hat = transpose(X_hat);


%%%% Calculating Symbol Errors for ZERO- FILTERING %%%%%%%%%%
number(a) = symerr(X_hat,x);

%%%%% Symbol Error Rate(SER)for ZERO- FILTERING %%%%%%%%%%%%%
SER(a) = number(a)/(N*L);

%%%%%%%%%%%%%%%%%%%%% Weiner Filtering %%%%%%%%%%%%%%%%%%%%%%%%
%%%%% RXY %%%%%
c_values= [1, 0.3,0.1,0.01];
RXY =zeros(N,1);
rx = zeros(N,1);
for k = 1:length(N)
        if k == 1
            rx(k) = 2*(A^2);
            RXY(k) = c_values(1)*rx(k);
            
        end
end

%%%%%%%% RY %%%%%%
Ry =zeros(N,N);
% Autocorrelation of noise signal Rv=1 if K=0, Rv=0 when k>1
rv =1;
for k = 1:N
   
         Ry(k,k) = c_values(1)*RXY(1)+rv;
         Ry(k,k+1) = c_values(2)*RXY(1);
         Ry(k+1,k) = c_values(2)*RXY(1);
         Ry(k,k+2) = c_values(3)*RXY(1);
         Ry(k+2,k) = c_values(3)* RXY(1);
         Ry(k,k+3) = c_values(4) *RXY(1);
         Ry(k+3,k) = c_values(4)*RXY(1);
end
% Removing 3 columns and 3 rows from Ry
Ry(N+1,:) =[];
Ry(N+2,:) =[];
Ry(N+1,:) =[];
Ry(:,N+1) =[];
Ry(:,N+2) = [];
Ry(:,N+1) = [];

%%%%%%% Weiner - Hopf Equations %%%%
W_w = inv(Ry)*RXY;


%%%%% Applying Weiner filter for Y generation  for WEINER-FILTER%%%%
Y_dash_W = conv(Y,W_w);

%%%% Finding Min A of Closest X to Y'(n) for WEINER-FILTER %%%%%%%%%
for o = 1:N*L
    for p = 1:length(QPSK_Set)
    
       
        all_Diff_w = Y_dash_W(o) - QPSK_Set(p);
        norm_diff_w = norm(all_Diff_w,2);
        norm_Store_w(p) = norm_diff_w;
    
     
    end
    [minvalue_w , minindx_w] = min(norm_Store_w);
     value_w(o) = minvalue_w;
     indx_w(o)= minindx_w;
    
end

%%%%%%%%%%%% X Hat for WEINER-FILTER %%%%%%%%%
X_hat_w = QPSK_Set(indx_w);
X_hat_w = transpose(X_hat_w);

%%%%%%% Calculating Symbol Errors for WEINER-FILTER %%%%%%%
number_w(a) = symerr(X_hat_w,x);

%%%%%%%%%%%%%% Symbol Error Rate(SER)for WEINER-FILTER %%%%%%%%%%%%
SER_W(a) = number_w(a)/(N*L);



end




SER
SER_W
close
figure
semilogy(SNR,SER,'-x','Linewidth',1);
hold on
semilogy(SNR,SER_W,'-o','Linewidth',1);
%axis([0 14 10^-5 0.5])
grid on
legend('Zero-Filtering','Weiner-Filtering','Location','southwest');
xlabel('SNR');
ylabel('SER');
title('Comparision of ZF and Wiener filter for QPSK in ISI channelling');
