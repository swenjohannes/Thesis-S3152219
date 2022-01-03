% Read vol data
filename = 'vol_surface_2012.csv';
vol_tbl = readtable(filename);
K = vol_tbl{1,2:end}/10;
Tvec = str2tenor(vol_tbl{2:end,1});
T = Tvec *[1,1/12,1/365,0,0,0]';
vol_dat = vol_tbl{2:end,2:end}/100;

figure(1)
[X, Y] = meshgrid(K, T);
%scatter3(X(:), Y(:), vol_dat(:));
surf(X, Y, vol_dat);
xlabel('Strike')
ylabel('Maturity')
zlabel('Volatility')

% Data
%S0 = 1.1325; %Spot
S0 = 1.1225; %Spot
r = 0;
q = 0;
v0 = 0.06^2;
kappa = 1.55;
theta = v0*2; 
eta = 0.17;
rho = -0.2;
x0 = [v0,kappa,theta,eta,rho];
H = 0.1;
pc = 1;
Xhs = zeros(length(T),5);
Xhr = zeros(length(T),5);

for j = 1:length(T)
    Tc = T(j);
    disp(Tc)
    vol = vol_dat(j,:);
    err_func = @(x) sum(sum((vol - get_iv(K,pc,Tc, S0,r,q,x)).^2));
    err_func_hr = @(x) sum(sum((vol - get_ivhr(K,pc,Tc, S0,r,q,x,H)).^2));

    lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
    ub = [1, 3, 10*v0, 1, 1];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    %get_iv(K,pc,Tc, S0,r,q,x0)
    options=optimset('Display','off');
    x_opt = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub,[],options);
    Xhs(j,:) = x_opt;
    x_opt = fmincon(err_func_hr,x0,[],[],[],[],lb,ub,[],options);
    get_iv(K,pc,Tc, S0,r,q,x0)
    Xhr(j,:) = x_opt;
end

%% 
figure(2)
%hold on
col = hsv(length(T));
for j = 1:length(T)
    plot(K,vol_dat(j,:),'o', K,get_iv(K,pc,T(j), S0,r,q,Xhs(j,:)),'Color',col(j,:))
    hold all
end
xlabel('Strike')
ylabel('Vol')
title('Classic Heston per maturity')
hold off

figure(3)
%hold on
col = hsv(length(T));
for j = 1:length(T)
    plot(K,vol_dat(j,:),'o', K,get_ivhr(K,pc,T(j), S0,r,q,Xhr(j,:),H),'Color',col(j,:))
    hold all
end
xlabel('Strike')
ylabel('Vol')
title('Rough Heston per maturity')
hold off
%% Glob
vol = vol_dat;
err_func = @(x) sum(sum((vol - get_iv(K,pc,T, S0,r,q,x)).^2));

lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
ub = [1, 3, 10*v0, 1, 1];
A = [];
b = [];
Aeq = [];
beq = [];
%get_iv(K,pc,Tc, S0,r,q,x0)
options = optimset('Display','off');
x_opt = fmincon(err_func,x0,A,b,Aeq,beq,lb,ub,[],options);

%% 
figure(4)
%hold on
col = hsv(length(T));
for j = 1:length(T)
    plot(K,vol_dat(j,:),'o', K,get_iv(K,pc,T(j), S0,r,q,x_opt),'Color',col(j,:))
    hold all
end
xlabel('Strike')
ylabel('Vol')
title('Classic Heston global')
hold off

%% Global fit rough Heston
vol = vol_dat;
err_func = @(x) sum(sum((vol - get_ivhr(K,pc,T, S0,r,q,x,H)).^2));

lb = [0.0001,0.0001, 0.0001, 0.0001, -1];
ub = [1, 3, 10*v0, 1, 1];

%get_iv(K,pc,Tc, S0,r,q,x0)
options = optimset('Display','off');
x_opt = fmincon(err_func,x0,[],[],[],[],lb,ub,[],options);

%% 
figure(5)
%hold on
col = hsv(length(T));
for j = 1:length(T)
    plot(K,vol_dat(j,:),'o', K,get_ivhr(K,pc,T(j), S0,r,q,x_opt,H),'Color',col(j,:))
    hold all
end
xlabel('Strike')
ylabel('Vol')
title('Rough Heston global')
hold off

Xhr_calib = [0.004114206022302   2.519859649229445   0.004080251303929   0.180672717676118  -0.405461005937409
   0.002818446719375   1.236083598865331   0.007671888654723   0.093208835579095  -0.234651427995674
   0.002547705430396   1.216664380498212   0.008804364864151   0.069540221489395  -0.118421784941760
   0.002951524361074   1.556507997361976   0.006452797457238   0.091232141287045  -0.137777051529585
   0.003234980868330   1.558342709634217   0.005688358069470   0.101744693315629  -0.164090630831127
   0.004292384866300   1.730552721994424   0.004554893595700   0.141471540323999  -0.262977845378220
   0.005056928361927   1.839207020487802   0.003896718864083   0.168134224808616  -0.262684094522906
   0.005912597433742   1.881195728930244   0.003698972225466   0.201716267051870  -0.277560216085186
   0.009430092068474   1.920051814163361   0.003721896191240   0.290001698124578  -0.275562913052522
   0.013398967139939   1.953003805422445   0.003781822384092   0.377412439147766  -0.270376022794821
   0.017983899318804   1.518266290913972   0.004613638834624   0.500493748900409  -0.260073626205014
   0.023740729077386   1.550774031825686   0.004727217246630   0.582108058261936  -0.276338982651167
   0.029803380371013   1.559959775505497   0.004847532975499   0.615780056529312  -0.299226842167389];

Xhs_calib = [0.003721931005252   2.849624282430895   0.016525924640131   0.984507534769602  -0.427549993262245
   0.002727669460471   1.544818921022308   0.023909650638108   0.411602192982312  -0.250020650898882
   0.002951111407622   1.302504732617253   0.014751765059156   0.209818211314443  -0.126816871286378
   0.002884436924011   1.274719192767200   0.013245565750325   0.187877991442363  -0.149184182200181
   0.002795938314238   1.250143405531782   0.011489322238670   0.172346804012891  -0.178016234304836
   0.003150875480933   1.368382557776410   0.007523776059115   0.167169691100250  -0.280684823961648
   0.003838450852050   1.660123573877532   0.005215933508212   0.181521867358402  -0.269842444624208
   0.004916868062789   1.845036197207846   0.004231327680092   0.191372033934980  -0.287021543825361
   0.011377272665497   2.052912153549729   0.003551971401226   0.265206117532209  -0.278527744676048
   0.019115192572559   1.973713495162116   0.003625811004416   0.347639885759077  -0.266727572763123
   0.032357325692205   1.516379720006511   0.004253241061877   0.459454643732245  -0.254732060056204
   0.051241183229493   1.551578699629970   0.004457862394576   0.567916992410410  -0.265806638955679
   0.076682826141243   1.591441644234797   0.004598993392917   0.614179168479151  -0.287462084448973];

