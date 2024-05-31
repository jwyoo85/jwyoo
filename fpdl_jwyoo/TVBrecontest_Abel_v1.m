clear all;
close all;

% Abel inversion for toroidal visible bremsstrahlung
% numerical abel inversion reference : http://zfn.mpdl.mpg.de/data/Reihe_A/46/ZNA-1991-46a-0639.pdf

%phantom_name = 'vb_KSTAR';
phantom_name = 'vb_JIPP_T-II';
%phantom_name = 'vb_EAST';
%phantom_name = 'vb_LHD';

load(['Phantom\',phantom_name,'.mat']);
if strcmp(phantom_name,'vb_JIPP_T-II')==1
    data = data1;
elseif strcmp(phantom_name,'vb_LHD')==1
    data = data2;
end

%% Abel inversion input
lch = 1;             % array number = 1
lsq = 1;             % 최소자승법 쓰려면 1 입력. 대수적으로 풀려면 1이 아닌값 입력
upf = 9;             % expansion 갯수. 일종의 최적화 파라미터 역할. 보통 4~30 사이에서 rms 가 최소값을 보여줌.
noise = 5;           % noise level (%)
N_test = 100;         % # of test for random noise

%% Reconstruction grid (length : mm unit)
R1 = 1300;           % radius of inner boundary (mm) (azimuthal sym. plasma)
R2 = 2300;           % radius of outer boundary
Rc = 1800;           % radius of center

angle = [0:0.01:2*pi]';
wall_in(:,1) = R1*cos(angle);   % inner boundary x
wall_in(:,2) = R1*sin(angle);   % inner boundary y
wall_out(:,1) = R2*cos(angle);  % outer boundary x
wall_out(:,2) = R2*sin(angle);  % outer boundary y
center(:,1) = Rc*cos(angle);  % plasma center x
center(:,2) = Rc*sin(angle);  % plasma center y

% line of sight of TVB
load('TVB_LoS\TVB_LOS.mat');
ch = 1:10;
badch = [];
%badch = [7,8];
if numel(badch)>0
    ch(badch) = [];
    los_rtheta(badch,:) = [];
    los_spot(badch,:) = [];
    los_xy(badch,:) = [];
end
los_rtheta(find(los_rtheta(:,2)>0),2) = los_rtheta(find(los_rtheta(:,2)>0),2)-180;

% los_rtheta : tangential position of each LoS (R, theta)
% los_xy1 : tangential position of each LoS (x, y)
los_xy1 = zeros(size(los_rtheta));
for i=1:size(los_xy1,1)
    los_xy1(i,1) = los_rtheta(i,1)*cos(los_rtheta(i,2)*pi/180);
    los_xy1(i,2) = los_rtheta(i,1)*sin(los_rtheta(i,2)*pi/180);
end

% vec : vector of Line of Sight
vec = zeros(size(los_xy1));
for i=1:size(vec,1)
    vec(i,1) = los_xy1(i,2);
    vec(i,2) = -los_xy1(i,1);
    vec(i,:) = vec(i,:)/norm(vec(i,:));
end

% los_xy : calibration spot position
% los_start/end : start & end position of Line of Sight
los_length = 4100;
los_start = los_xy - 0.5*los_length*vec;
los_end = los_xy + 0.5*los_length*vec;

%% Reconstruction grid (mm)
NR = 10;
%R_min = 1800;       R_max = 2300;
R_min = min(los_rtheta(:,1));  R_max = max(los_rtheta(:,1));
R = linspace(R_min,R_max,NR)+30;
dR = R(2)-R(1);

% load('Phantom\R_ts.mat');
% R = R_ts*1000;
% NR = length(R);
% R_min = min(R);
% R_max = max(R);

% R = transpose(flipud(los_rtheta(:,1)));
% NR = length(R); R_min = min(R); R_max = max(R);

rho = (R-1800)/500;
phantom = interp1(data(:,1),data(:,2),rho);

n_in = find(rho<data(1,1));
n_out = find(rho>data(end,1));
if numel(n_in)>0
    phantom(n_in) = phantom(max(n_in)+1);
end
if numel(n_out)>0
    phantom(n_out) = phantom(min(n_out)-1);
end

%% Plot of reconstruction grid & TVB LOS
xmin = -2300;
xmax = 2300;
ymin = -2300;
ymax = 2300;

figure; hold on;
plot(wall_in(:,1),wall_in(:,2),'-k');
plot(center(:,1),center(:,2),'--k');
plot(wall_out(:,1),wall_out(:,2),'-k');
for i=1:size(vec,1)
    x1 = [los_start(i,1), los_end(i,1)];
    y1 = [los_start(i,2), los_end(i,2)];
    plot(x1, y1, '-b'); hold on;
end
for i=1:length(R)
    x1 = R(i)*cos(angle);   % inner wall x
    y1 = R(i)*sin(angle);   % inner wall y
    plot(x1, y1, '-r'); hold on;
end
%scatter(los_xy1(:,1),los_xy1(:,2),'b');
%scatter(los_xy(:,1),los_xy(:,2),'r');
xlim([xmin xmax]);  ylim([ymin ymax]);
axis image; grid on;

%% calculation of cosine expansion
% Details: The unknown distribution f(r) is expanded as
%           f(r) = sum_{n=lof}^{upf} (A_n * f_n(r))                    (1)
% where the lower frequency is set to 1 and the upper frequency upf is
% important for noise-filtering. f_n(r) is a set of cos-functions:
%           f_n(r) = 1 - (-1)^n*cos(n*pi*r/R)  and f_0(r) = 1          (2)
% For the Abel inversion, the integrals h_n have to be calculated
%           h_n(d) = int_y^R f_n(r) * r / sqrt(r^2-d^2) dr             (3)

d1 = los_rtheta(:,1);   %
RM = R_max;     % reconstruction boundary
s = R;    % reconstruction grid
fn = zeros(length(s), upf+1);
hn = zeros(length(d1), upf+1);

% for n=0 case
fn(:,1) = 1;
for i=1:length(d1)
    d = d1(i);
    n0 = @(t) ones(size(t));
    hn(i,1) = integral(n0,0,sqrt(RM^2-d^2));
end

%for n>0 case
for n=1:upf
    for i=1:length(s)
        fn(i,n+1) = (1-(-1)^n*cos(n*pi*s(i)/RM)); % f(r) 함수의 코사인 익스팬션
    end
    for i=1:length(d1)
        d=d1(i);
        nn = @(t) (1-(-1)^n*cos(n*pi*sqrt(t.^2+d.^2)/RM)).*ones(size(t));
        hn(i,n+1) = integral(nn,0,sqrt(RM^2-d^2)); % h(r) 함수의 코사인 익스팬션
    end
end

%% line-integration
RM = R_max;     % reconstruction boundary
h_ph0 = zeros(size(d1));     % 선적분 신호

for i=1:length(h_ph0)
    n_int1 = find(R>d1(i),1,'first');
    n_int2 = find(R<RM,1,'last');
    phantom_temp = phantom(n_int1:n_int2);
    R_temp = R(n_int1:n_int2);
    h_ph0(i) = 2*sum(phantom_temp.*R_temp./sqrt(R_temp.^2-d1(i)^2)*dR);
end

h_ph = zeros(size(d1,1),N_test);
for i=1:N_test
    h_ph(:,i) = h_ph0 + max(h_ph0)*(rand(size(h_ph0))-0.5)*noise/100;
end

%% solve equation system {H = 2*sum_n (A_n*hn)}, or say {B = A*L} for the amplitude A
% 최소자승법: (sum_{k=1}^{N} (H(y_k) - h(y_k))^2 == min) has to be solved
% H(y): fn 함수의 적분꼴(numerical)
% h(y): f 함수 적분꼴(analytic)
% h_ph: phantom의 선적분 신호

A = zeros(upf+1,N_test);
for k=1:N_test
    if lsq == 1
        x0 = 1*ones(upf+1,1);    % guess some initial values for optimisation
        optfun = @(x,hn) 2*hn*x; % function to be optimized
        A(:,k) = lsqcurvefit(optfun,x0,hn,h_ph(:,k));
    else
        B = zeros(1,upf+1);
        L = zeros(upf+1, upf+1);
        for i=1:upf+1
            for j=1:upf+1
                L(j,i) = 2.*sum(hn(:,i).*hn(:,j));
            end
            B(i) = sum(h_ph(:,k).*hn(:,i));
        end
        A(:,k) = B/L;
    end
end

%% local emissivity profile 구하기
f_rec = zeros(length(s),N_test); % 재구성된 프로파일
for k=1:N_test
    f_rec(:,k) = f_rec(:,k) + A(1,k)*1;     % special case for n=0 where f_0(r) = 1
    % case for n>0
    for i=2:upf+1
        f_rec(:,k) = f_rec(:,k) + A(i,k).*fn(:,i);
    end
end

e_rms = zeros(1,N_test);
for k=1:N_test
    e_rms(k) = norm(phantom-f_rec(:,k))/norm(phantom)*100;
end
e_avg = mean(e_rms);

figure; hold on;
plot(R,phantom,'-b',R,f_rec(:,1),'-r');
xlabel('R (mm)');  ylabel('counts');
legend('phantom','recon.');