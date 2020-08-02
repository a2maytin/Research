function sigma_a = NormMSDInterceptErrorNW(x, P, N) 

sum0 = 0;
sum1 = 0;

if (P < 2 || P >= N)
      sigma_a =  0;
end
%for (i=1; i<= P; i++)
for i = 1:P
    fi = ff(i, N, x);
    sum0 = sum0 + ((2*P+1)/3 - i)^2/fi; 
    psum_i = 0;
    %for (j=1; j<i; j++)
    for j = 1:(i-1)
        gij = g(i,j,N,x);
        psum_i = psum_i + ((2*P+1)/3 - j)*gij;
    end
    sum1 = sum1 + psum_i*((2*P+1.)/3. - i);
end
temp = (sum0 + 2*sum1);
sigma_a = sqrt(temp)*6/P/(P-1)/x; 



function f = ff(n, N, x) 

K = N - n; 
if (n <= K)
    f = fminus(n,N,x);
else
    f = fplus(n,N,x);
end   
end

function fmin =  fminus(n, N, x)
K = N - n;
f = n*(4*n^2*K + 2*K - n^3 + n)/6/K^2 + (2*n*x + (1+(1 - ...
    n/K)/2)*x^2)/K;
fmin = 1/f;
end

function fplu =  fplus(n, N, x)
K = N - n;
f = (6*n^2*K - 4*n*K^2 + 4*n + K^3 - K)/6/K + (2*n*x + x^2)/K; 
fplu = 1/f;
end

function g = g(m, n, N, x) 
K = N-n;
P = N-m;
if (m <= n) 
      temporary = m;
      m = n;
      n = temporary;
end
if (m+n > N)
    g = (N^2*(N - 4*n - 3*m) + N*(6*n^2 + 8*m*n + 3*m^2 - 1) - ...
    6*n^3 - 4*m^2*n + 4*n - m^3 + m)/6/K ...
    + (2*n*x + x^2/2)/K;
else
    g = n*(-N*(2*n^2 - 6*m*n  -2) ...
    - n^3 + 2*m*n^2 - 6*m^2*n + n - 2*m)/6/K/P + ...
    (2*n*x +  (1 - n/2/P)*x^2/2)/K;   
end
end

end



