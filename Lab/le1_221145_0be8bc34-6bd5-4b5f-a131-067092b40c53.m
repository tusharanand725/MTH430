T = 1;

M = 10;
for m = 1:M
    N = 2^m;
    h = T/N;

    y = zeros(1,N+1);  % this array holds the exact solution
    yE = zeros(1,N+1); % holds the numerical solution obtained using Euler's method
    yB = zeros(1,N+1); % holds the numerical solution obtained using Backward Euler
    yT = zeros(1,N+1); % holds the numerical solution obtained using Trapezoidal method
   y(1)=2;
    yE(1)=2;
    yT(1)=2;
    yB(1)=2;
    t=linspace(0,1,N+1);
    %exact solution
    for i=1:N
        y(i+1)= exp(-50*t(i)) + exp(2*sin(20*t(i)));
    end
    for n= 1:N
        %fuctions for trapezoidal method
        f1 = (-50*(yT(n)-exp(2* sin(20*t(n)) ) ) + 40* cos(20*(t(n))) * exp(2*sin(20*t(n))) );
        f2 = ((50*exp(2* sin(20*t(n+1)) ) ) + 40* cos(20*(t(n+1))) * exp(2*sin(20*t(n+1))) );
        %euler method 
    yE(n+1)=yE(n) + h*(-50*(yE(n)-exp(2* sin(20*t(n)) ) ) + 40* cos(20*(t(n))) * exp(2*sin(20*t(n))) );
    %backward euler method
    yB(n+1)=(yB(n) + h*((50* exp(2* sin(20*t(n+1)) ) ) + 40* cos(20*(t(n+1))) * exp(2*sin(20*t(n+1))) ))/(1+50*h);
    %trapezoidal method
    yT(n+1)=(yT(n) + (h/2)*(f1 + f2 ))/(1+25*h);
    end
    % Do not change the template below this line
    %%

    eE(m) = max(abs(yE-y))/max(abs(y));
    eB(m) = max(abs(yB-y))/max(abs(y));
    eT(m) = max(abs(yT-y))/max(abs(y));
end

% Convergence report
for m = 2:M
    fprintf('%4d \t %.2e \t %.1e \t %.2e \t %.1e \t %.2e \t %.1e \n', ...
        2^m, eE(m), eE(m-1)/eE(m), eB(m), eB(m-1)/eB(m), ...
        eT(m), eT(m-1)/eT(m));
end

plot([0:N]*h,y,'LineWidth',2);
hold on;
plot([0:N]*h,yE,'--','LineWidth',2);
plot([0:N]*h,yB,'-.','LineWidth',2);
plot([0:N]*h,yT,':','LineWidth',2);
legend('y','yE','yBE','yT');
