%% Data
X = X;      % covariate vector
d_data = Y; % response vector
test_X = Test_X(:,1);
test_Y = Test_Y(:,1);

wd_data = sort(d_data);
sort_data = sortrows([X';d_data']');

%% Parameter
alpha = 0.5;     % coverage level
b = 20;          % number of polynomial pieces
start_value = 0; % covariate left endpoint
end_value = 10;  % covariate right endpoint
d = 2;           % degree of each polynomial
rho = 1;         % smoothness of spline
lambda = 0.001;  % smoothing parameter
iter = 1000;     % number of  ADMM iterations

h_x = (end_value - start_value) / b;
spline_knots = start_value:h_x:end_value;
n_s = length(X);
alpha_number = length(alpha);
L = alpha_number * 2;
[~, ~, bandwidth] = ksdensity(X, 'Bandwidth', 'plug-in');


%% Generate matriics
iiii=0;
kkkk=ones(1,b+1);
for kk=1:b
    for i = 1:n_s
        if sort_data(i,1) <= kk*1 && sort_data(i,1)>(kk-1)*1
            iiii = iiii+1;
        end
        kkkk(kk+1) = iiii;
    end    
end

%%%Weight Vector: w%%%
w = zeros(n_s,1);
for i =1:n_s
    w(i,1)=nthroot(ksdensity(sort_data(:,1),sort_data(i,1),'Bandwidth', 'plug-in'),5);
end

%%%Quantile Level Vector: p%%%
p_U = zeros(n_s,alpha_number);
p_L = zeros(n_s,alpha_number);
SORT = sortrows(sort_data,2);
for i = 1:n_s
    [p_L(i,:),p_U(i,:)] = cmie(SORT,sort_data(i,1),alpha,bandwidth);
end

%%%Observations Matrix: A%%%
iii = ones(L*n_s*(d+1),1);
jjj = ones(L*n_s*(d+1),1);
A_s = zeros(L*n_s*(d+1),1);
ii = 1;
for kk=1:L
    bk = 1;
    for i = 1:n_s
        if spline_knots(bk+1)< sort_data(i,1)
            bk = bk+1;
        end
        tau = (sort_data(i,1)-spline_knots(bk))/h_x;
        for dd = 0:d 
            iii(ii+dd) = (kk-1)*n_s+i;
            jjj(ii+dd) = (kk-1)*b*(d+1) + (bk-1)*(d+1)+dd+1;
            A_s(ii+dd) = tau^(d-dd);
        end
        ii = ii + d+1;
    end
end
A = sparse(iii,jjj,A_s,L*n_s,L*b*(d+1));

%%%Regularzation Term Matrix: Q%%%
Q_small=zeros(d-1,d-1);
for i=0:d-2
    for j=0:d-2
        Q_small(i+1,j+1)=(d-i)*(d-i-1)*(d-j)*(d-j-1)/(2*d-i-j-3);
    end
end
Q_small=(Q_small+Q_small')/2;

iii=zeros(L*b*((d-1)^2),1);
jjj=zeros(L*b*((d-1)^2),1);
Q_s=zeros(L*b*((d-1)^2),1);
for kk=1:L
    for ii=1:b
        for i=0:d-2
            for j=0:d-2
                iii((kk-1)*b*((d-1)^2) + (ii-1)*((d-1)^2)+i*(d-1)+j+1)=(kk-1)*b*(d+1) + (ii-1)*(d+1)+i+1;
                jjj((kk-1)*b*((d-1)^2) + (ii-1)*((d-1)^2)+i*(d-1)+j+1)=(kk-1)*b*(d+1) + (ii-1)*(d+1)+j+1;
                Q_s((kk-1)*b*((d-1)^2) + (ii-1)*((d-1)^2)+i*(d-1)+j+1)=(1/(h_x^3))*Q_small(i+1,j+1);
            end
        end
    end
end
Q=sparse(iii,jjj,Q_s,L*b*(d+1),L*b*(d+1));
Q=(Q+Q')/2;

%%%Differentability Matrix: H%%%
iii=zeros(L*(b-1)*(2*(d+2)-rho)*(rho+1)/2,1);
jjj=zeros(L*(b-1)*(2*(d+2)-rho)*(rho+1)/2,1);
H_s=zeros(L*(b-1)*(2*(d+2)-rho)*(rho+1)/2,1);
for kk=1:L
    for ii=1:b-1
        for l=0:rho
            for k=l:d
                iii((kk-1)*(b-1)*(2*(d+2)-rho)*(rho+1)/2 + (ii-1)*(2*(d+2)-rho)*(rho+1)/2+(2*(d+2)-l+1)*l/2+k-l+1)=(kk-1)*(rho+1)*(b-1) + (ii-1)*(rho+1)+l+1;
                jjj((kk-1)*(b-1)*(2*(d+2)-rho)*(rho+1)/2 + (ii-1)*(2*(d+2)-rho)*(rho+1)/2+(2*(d+2)-l+1)*l/2+k-l+1)=(kk-1)*b*(d+1) + (ii-1)*(d+1)+d-k+1;
                H_s((kk-1)*(b-1)*(2*(d+2)-rho)*(rho+1)/2 + (ii-1)*(2*(d+2)-rho)*(rho+1)/2+(2*(d+2)-l+1)*l/2+k-l+1)=(1/(h_x^l))*factorial(k)/factorial(k-l);
            end
            iii((kk-1)*(b-1)*(2*(d+2)-rho)*(rho+1)/2 + (ii-1)*(2*(d+2)-rho)*(rho+1)/2+(2*(d+2)-l+1)*l/2+d-l+2)=(kk-1)*(rho+1)*(b-1) + (ii-1)*(rho+1)+l+1;
            jjj((kk-1)*(b-1)*(2*(d+2)-rho)*(rho+1)/2 + (ii-1)*(2*(d+2)-rho)*(rho+1)/2+(2*(d+2)-l+1)*l/2+d-l+2)=(kk-1)*b*(d+1) + ii*(d+1)+d-l+1;
            H_s((kk-1)*(b-1)*(2*(d+2)-rho)*(rho+1)/2 + (ii-1)*(2*(d+2)-rho)*(rho+1)/2+(2*(d+2)-l+1)*l/2+d-l+2)=-(1/(h_x^l))*factorial(l);
        end
    end
end
H=sparse(iii,jjj,H_s,L*(rho+1)*(b-1),L*b*(d+1));

%%%Noncrossing Matrix: G%%%
iii=ones((L-1)*(b*(d+1)*(d+2)),1);
jjj=ones((L-1)*(b*(d+1)*(d+2)),1);
G_s=zeros((L-1)*(b*(d+1)*(d+2)),1);
for kk=1:(L-1)
    for ii=1:b
        for l=0:d
            for k=0:l
                iii((kk-1)*(b*(d+1)*(d+2)) + (ii-1)*(d+1)*(d+2)+l*(l+1)+2*k+1)=(kk-1)*b*d + (ii-1)*d+l+1;
                jjj((kk-1)*(b*(d+1)*(d+2)) + (ii-1)*(d+1)*(d+2)+l*(l+1)+2*k+1)=(kk-1)*b*(d+1) + (ii-1)*(d+1)+d+1-k;
                G_s((kk-1)*(b*(d+1)*(d+2)) + (ii-1)*(d+1)*(d+2)+l*(l+1)+2*k+1)=-factorial(d-k)/(factorial(l-k)*factorial(d-l));

                iii((kk-1)*(b*(d+1)*(d+2)) + (ii-1)*(d+1)*(d+2)+l*(l+1)+2*k+2)=(kk-1)*b*d + (ii-1)*d+l+1;
                jjj((kk-1)*(b*(d+1)*(d+2)) + (ii-1)*(d+1)*(d+2)+l*(l+1)+2*k+2)=kk*b*(d+1) + (ii-1)*(d+1)+d+1-k;
                G_s((kk-1)*(b*(d+1)*(d+2)) + (ii-1)*(d+1)*(d+2)+l*(l+1)+2*k+2)=factorial(d-k)/(factorial(l-k)*factorial(d-l));
            end
        end
    end
end
G=sparse(iii,jjj,G_s,(L-1)*b*(d+1),L*b*(d+1));

%%% Solved by ADMM %%%
gamma = 1;
[LLL1,UUU1] = lu(2*lambda*Q + (A'*A+ G'*G)/gamma);
[LLL2,UUU2] = lu( H*(UUU1\(LLL1\(H'))));

c_k = rand(L*b*(d+1), 1);
z1_k = A*c_k;
z2_k = G*c_k;

u1_k = z1_k;
u2_k = z2_k;
dif1 = 1;
dif2 = 1;
dif3 = 1;
iteration = 1;
while dif1 > 10^(-5) || dif2 > 10^(-5) || dif3 > 10^(-5)
    iteration = iteration + 1;
    aaa = (A'*(z1_k-u1_k)+G'*(z2_k-u2_k))/gamma;
    c_k1 = UUU1\(LLL1\(aaa - H'*(UUU2\(LLL2\(H*(UUU1\(LLL1\aaa)))))));
    z1_k1 = A*c_k1 + u1_k;

    for II = 1:alpha_number
        for I = 1:n_s
            z1_k1(n_s * (II-1) + I) = per_prox(z1_k1(n_s * (II-1) + I) ,sort_data(I,2) ,p_L(I,II),gamma*w(I));
        end
    end
    for JJ = 1:alpha_number
        for J = 1:n_s
            z1_k1(n_s * (alpha_number + JJ - 1) + J) = per_prox(z1_k1(n_s * (alpha_number + JJ - 1) + J) ,sort_data(J,2) ,p_U(J,alpha_number - JJ + 1),gamma*w(J));
        end
    end

    z2_k1 = G*c_k1 + u2_k;
    z2_k1(z2_k1<0)=0;

    u1_k1 = u1_k + A*c_k1 -  z1_k1;
    u2_k1 = u2_k + G*c_k1 -  z2_k1;

    dif1 = norm(c_k1-c_k);
    dif2 = norm([z1_k1;z2_k1]-[z1_k;z2_k]);
    dif3 = norm([u1_k1;u2_k1]-[u1_k;u2_k]);
    c_k = c_k1;
    z1_k = z1_k1;
    z2_k = z2_k1;
    u1_k = u1_k1;
    u2_k = u2_k1;
    if iteration>iter
        break
    end 
end


%% Output
c_star = c_k;
fineness = 100;
xxxx = spline_knots(1):h_x/fineness:spline_knots(end);
zzzz = zeros(fineness*b+1,L);

for iii=1:L
    for ii=1:b
        for i=0:(fineness-1)
            tau=(xxxx((ii-1)*fineness+i+1)-spline_knots(ii))/h_x;
            for k=0:d
                zzzz((ii-1)*fineness+i+1,iii) = zzzz((ii-1)*fineness+i+1,iii) + c_star((iii-1)*(b*(d+1))+(ii-1)*(d+1)+d-k+1)*tau^k;
            end
        end
    end
    for k=0:d
        zzzz(fineness*b+1,iii) = zzzz(fineness*b+1,iii) + c_star((iii-1)*(b*(d+1))+(b-1)*(d+1)+d-k+1);
    end
end

middle_curve = (zzzz(:,L/2) + zzzz(:,L/2+1))/ 2;

figure(1)
hold on

for i = 1:L/2
    fill([xxxx,fliplr(xxxx)],[zzzz(:,i)',fliplr(zzzz(:,L-i+1)')],[0.8500, 0.3250, 0.0980]); %estimated MI
end
plot(xxxx,zzzz(:,1),'k','LineWidth',5);                                                     % upper bound
plot(xxxx,zzzz(:,2),'k','LineWidth',5);                                                     % lower bound
plot(xxxx,middle_curve,':','Color','black','LineWidth',5);                                  % middle curve
hold on

scatter(X',d_data',12,'blue','filled');
xlim([start_value end_value]);
set(gca, 'FontSize', 20, 'LineWidth', 3);
xlabel('X');
ylabel('Y');

save('mir_coefficient.mat','c_star');

%% mCWC
n = length(test_X);
lower = zeros(n, 1);
upper = zeros(n, 1);

nCovered = 0;
sumWidth = 0;

c_star_m1 = c_star(1 : b*(d+1));
c_star_m2 = c_star(b*(d+1)+1 : 2*b*(d+1));

for i = 1:n
    ii = 0;
    for j = start_value : h_x : (end_value - h_x)
        ii = ii + 1;
        if (j <= test_X(i)) && (test_X(i) < j + h_x)
            tau = (test_X(i) - j) / h_x;
            break
        end
    end

    for k = 0:d
        idx = (ii - 1) * (d + 1) + (d - k + 1);
        lower(i) = lower(i) + c_star_m1(idx) * tau^k;
        upper(i) = upper(i) + c_star_m2(idx) * tau^k;
    end

    if (test_Y(i) <= upper(i)) && (test_Y(i) >= lower(i))
        nCovered = nCovered + 1;
    end

    sumWidth = sumWidth + (upper(i) - lower(i));
end

PICP = nCovered / n;
MPIW = sumWidth / n;

yRange = max(test_Y) - min(test_Y);
NMPIW = MPIW / yRange;

if PICP < alpha
    mCWC_value = NMPIW * exp(-eta * (PICP - alpha));
else
    mCWC_value = NMPIW;
end

fprintf('mCWC = %.6f\n', mCWC_value);