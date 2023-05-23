% my student ID: 4834186, thus:
a = 4; b = 3; c = 6;

A = [a - b, 0.5 - c; 0, 1];
B = [0; 1];

%% question 1 - 1
p1 = [-1 + 2j, -1 - 2j];
p2 = [-1, -2];
p3 = [1, 2];
% format rat;
K1 = place(A, B, p1);
K2 = place(A, B, p2);
K3 = place(A, B, p3);

%% question 1 - 2
syms h s tau

Fx_h = expm(A * h);
Fu_h = int(expm(A*s)*B, s, [h-tau, h]); %zeros(2, 1);
G1 = int(expm(A*s)*B, s, [0 h-tau]);

F = [Fx_h, Fu_h; zeros(1,3)];
G = [G1; eye(1)];
G1

tau = 0;
F = subs(F);
G = subs(G);

hs = 0.01:0.01:0.6;
rho = zeros(3, length(hs));
cnt = 1;
for st=hs
    Fn = subs(F, h, st);
    Gn = subs(G, h, st);
    rho(1, cnt) = vpa(max(abs(eig(Fn - Gn * [K1, 0]))));
    rho(2, cnt) = vpa(max(abs(eig(Fn - Gn * [K2, 0]))));
    rho(3, cnt) = vpa(max(abs(eig(Fn - Gn * [K3, 0]))));
    cnt = cnt + 1;
end
figure(1);
plot(hs, rho(1,:), hs, rho(2,:), hs, rho(3,:), hs, ones(length(hs)), "--", [0.42, 0.575], [1, 1], "ro")
xlabel("sampling time h[s]")
ylabel("max. eigenvalue distance from origin")
legend("sys 1", "sys 2", "sys 3")
saveas(gcf,'figures4report/q1-2.png')



% sim_time = 3;
% T = sim_time/h;
% x0 = [0.1; 0];
% u = zeros(1, T);
% xe = [[x0; u(1)], zeros(3, T-1)];
% 
% 
% 
% 
% for k=1:1:T
% %     xe(:,k + 1) = F * xe(:,k) + G * u(k);
%     xe(:,k + 1) = (F - G * [K1, 0]) * xe(:,k);
% end
% 
% plot(1:T+1, xe(1,:), 1:T+1, xe(2,:))
% legend("x1", "x2")
% ylabel("state value")
% xlabel("time step")


%% question 2
syms h s tau

Fx_h = expm(A * h);
Fu_h = int(expm(A*s)*B, s, [h-tau, h]);
G1 = int(expm(A*s)*B, s, [0 h-tau]);

F = [Fx_h, Fu_h; zeros(1,3)];
G = [G1; eye(1)];

%% 


hs = 0.01:0.03:0.6;
taus = 0:0.03:0.55;
rho = zeros(length(taus), length(hs));
rho2 = zeros(length(taus), length(hs));

cnt_tau = 1;
for tau_i=taus
    cnt_h = 1;
    for st=hs
        Fn = subs(F, [h, tau], [st, tau_i]);
        Gn = subs(G, [h, tau], [st, tau_i]);
        rho(cnt_tau, cnt_h) = real(vpa(max(abs(eig(Fn - Gn * [K1, 0])))));
        rho2(cnt_tau, cnt_h) = real(vpa(max(abs(eig(Fn - Gn * [K2, 0])))));
        cnt_h = cnt_h + 1;
    end
    cnt_tau = cnt_tau + 1;
end

figure(2)
redChannel = rho >= 1;greenChannel = rho < 1;blueChannel=zeros(size(greenChannel));
C = double(cat(3, redChannel, greenChannel, blueChannel));
surf(hs, taus,rho,C)
xlabel("sampling time h[s]")
ylabel("delay time \tau[s]")
zlabel("\rho")
saveas(gcf,'figures4report/q2-1a.png')

figure(3)
redChannel = rho2 >= 1;greenChannel = rho2 < 1;blueChannel=zeros(size(greenChannel));
C = double(cat(3, redChannel, greenChannel, blueChannel));
surf(hs, taus,rho2,C)
xlabel("sampling time h[s]")
ylabel("delay time \tau[s]")
zlabel("\rho")
saveas(gcf,'figures4report/q2-1b.png')

%% Question 2-2 designing the dynamic controller - trial and error

h1 = 0.5; %sampling time for K1, stable with zero Tau
% h2 = ;

taus = 0:0.025:0.4;
rho1 = zeros(1, length(taus));
rho2 = zeros(1, length(taus));
rho3 = zeros(1, length(taus));

% rho2 = zeros(length(taus), length(hs));

cnt_tau = 1;
for tau_i=taus
    Fn = subs(F, [h, tau], [h1, tau_i]);
    Gn = subs(G, [h, tau], [h1, tau_i]);
    rho1(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * [K1, 0])))));
    rho2(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * [K1, 0.15])))));
    cnt_tau = cnt_tau + 1;
end
figure(4)
set(gca,'ColorOrderIndex',1)
plot(taus, rho1, taus, rho2)
xlabel("delay time \tau[s]")
ylabel("\rho")
legend("K_u = 0", "K_u = nonzero")
saveas(gcf,'figures4report/q3-1-stupid.png')


%% Q2-2 smarter approach - sys1
h1 = 0.5;
tau1 = 0.15;
p1 = [0.3, 0.2, 0.1];
p2 = [0.5, 0.6, 0.7];


Fnew = double(subs(F, [h, tau], [h1, tau1]));
Gnew = double(subs(G, [h, tau], [h1, tau1]));
Knew1 = vpa(place(Fnew, Gnew, p1));
Knew2 = vpa(place(Fnew, Gnew, p2));


taus = 0:0.025:0.4;
rho0 = zeros(1, length(taus)); % benchmark
rho1 = zeros(1, length(taus));
rho2 = zeros(1, length(taus));

% rho2 = zeros(length(taus), length(hs));

cnt_tau = 1;
for tau_i=taus
    Fn = subs(F, [h, tau], [h1, tau_i]);
    Gn = subs(G, [h, tau], [h1, tau_i]);
    rho0(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * [K1, 0])))));
    rho1(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * Knew1)))));
    rho2(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * Knew2)))));
    cnt_tau = cnt_tau + 1;
end
figure(5)
set(gca,'ColorOrderIndex',1)
plot(taus, rho0, taus, rho1, taus, rho2)
xlabel("delay time \tau[s]")
ylabel("\rho")
legend("K_u = 0", "K_u via p1", "K_u via p2")
saveas(gcf,'figures4report/q3-1-place-sys1.png')

%% Q2-2 smarter approach - sys2
h1 = 0.3;
tau1 = 0.2;
p1 = [0.3, 0.2, 0.1];
p2 = [0.5, 0.6, 0.7];


Fnew = double(subs(F, [h, tau], [h1, tau1]));
Gnew = double(subs(G, [h, tau], [h1, tau1]));
Knew1 = vpa(place(Fnew, Gnew, p1));
Knew2 = vpa(place(Fnew, Gnew, p2));


taus = 0:0.025:0.4;
rho0 = zeros(1, length(taus)); % benchmark
rho1 = zeros(1, length(taus));
rho2 = zeros(1, length(taus));

% rho2 = zeros(length(taus), length(hs));

cnt_tau = 1;
for tau_i=taus
    Fn = subs(F, [h, tau], [h1, tau_i]);
    Gn = subs(G, [h, tau], [h1, tau_i]);
    rho0(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * [K2, 0])))));
    rho1(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * Knew1)))));
    rho2(cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * Knew2)))));
    cnt_tau = cnt_tau + 1;
end
figure(6)
set(gca,'ColorOrderIndex',1)
plot(taus, rho0, taus, rho1, taus, rho2)
xlabel("delay time \tau[s]")
ylabel("\rho")
legend("K_u = 0", "K_u via p1", "K_u via p2")
saveas(gcf,'figures4report/q3-1-place-sys2.png')





%% Question 3 - 1
syms h s tau

% int 2 is
% int(expm(A*s)*B, s, [h-tau, h]);
Fu1 = int(expm(A*(h-s))*B, s, [-h+tau,tau]);
Fu2 = int(expm(A*(h-s))*B,s,[0, tau-h]);


Fx_h = expm(A * h);
Fu_h = int(expm(A*s)*B, s, [h-tau, h]);
G1 = int(expm(A*s)*B, s, [0 h-tau]);

F = [Fx_h, Fu_h; zeros(1,3)];
G = [G1; eye(1)];

hs = 0.1:0.025:0.5;
taus = 0:0.025:0.5;
rho = zeros(length(taus), length(hs));

cnt_tau = 1;
for tau_i=taus
    cnt_h = 1;
    tau_i
    for h_i=hs
        h_i
        if tau_i < h_i
            Fn = [Fx_h, Fu_h, zeros(2,1); 0,0,0,0;0,0,1,0];
            Gn = [G1; 1; 0];
        else
            Fn = [Fx_h, Fu1, Fu2; 0,0,0,0;0,0,1,0];
            Gn = [0;0; 1; 0];
        end
        Fn = subs(Fn, [h, tau], [h_i, tau_i]);
        Gn = subs(Gn, [h, tau], [h_i, tau_i]);
        rho(cnt_tau, cnt_h) = real(vpa(max(abs(eig(Fn - Gn * [K1, 0, 0])))));
        cnt_h = cnt_h + 1;
    end
    cnt_tau = cnt_tau + 1;
end

%% Q3 plotting stability

figure(7)
redChannel = rho >= 1;greenChannel = rho < 1;blueChannel=zeros(size(greenChannel));
C = double(cat(3, redChannel, greenChannel, blueChannel));
surf(hs, taus,rho,C)
hold on
w=gobjects(2,1);
w(1)=plot3(nan,nan,nan,'r*');
w(2)=plot3(nan,nan,nan,'g*');
legend(w, ["unstable" "stable"], 'Location','northeast')
xlabel("sampling time h[s]")
ylabel("delay time \tau[s]")
zlabel("\rho")
saveas(gcf,'figures4report/q33-stability.png')

%% Q3 controller design

h_g = 0.2;
tau_g = 0.29;
% p_g = [0.5+0.5j, 0.5-0.5j, 0.6+0.5j, 0.6-0.5j]; best so far
p_g = [0.6+0.51j, 0.6-0.51j, 0.6+0.5j, 0.6-0.5j];


Fnn = [Fx_h, Fu1, Fu2; 0,0,0,0;0,0,1,0];
Gnn = [0;0; 1; 0];

Fnn = subs(Fnn, [h, tau], [h_g, tau_g]);
Gnn = subs(Gnn, [h, tau], [h_g, tau_g]);

K_g = place(double(Fnn), double(Gnn), p_g);

taus = 0:0.05:0.3;

rho = zeros(1, length(taus));
rho_g = zeros(1, length(taus));


cnt_tau = 1;
for tau_i=taus
    if tau_i < h_g
        Fn = [Fx_h, Fu_h, zeros(2,1); 0,0,0,0;0,0,1,0];
        Gn = [G1; 1; 0];
    else
        Fn = [Fx_h, Fu1, Fu2; 0,0,0,0;0,0,1,0];
        Gn = [0;0; 1; 0];
    end
    Fn = subs(Fn, [h, tau], [h_g, tau_i]);
    Gn = subs(Gn, [h, tau], [h_g, tau_i]);
    rho(1, cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * [K1, 0, 0])))));
    rho_g(1, cnt_tau) = real(vpa(max(abs(eig(Fn - Gn * K_g)))));
    cnt_h = cnt_h + 1;
    cnt_tau = cnt_tau + 1;
end

figure(8)
plot(taus, rho, taus, rho_g)
legend("static controller", "dynamic controller")
xlabel("delay time \tau [s]")
ylabel("\rho")
saveas(gcf,'figures4report/q3-2.png')


%% Question 4 controller design
syms h s tau
K = [K1, 1];
Fx_h = expm(A * h);
Fu_h = int(expm(A*s)*B, s, [h-tau, h]);
G1 = int(expm(A*s)*B, s, [0 h-tau]);

F = [Fx_h, Fu_h; zeros(1,3)];
G = [G1; eye(1)];
ks = 0:0.05:1;
hns = 0:0.025:0.7;


rho = zeros(length(ks), length(hns));

cnt_k = 1;
for k_u=ks
    cnt_h = 1;
    for hn=hns
        h_i = 1.5 * hn + 0.000001;
        tau_i = 0.5* h_i;
        Fn = subs(F, [h, tau], [h_i, tau_i]);
        Gn = subs(G, [h, tau], [h_i, tau_i]);
        Kn = [K1, k_u];
        rho(cnt_k, cnt_h) = real(vpa(max(abs(eig(Fn - Gn * Kn)))));
        cnt_h = cnt_h + 1;
    end
    cnt_k = cnt_k +1;
end

%% plotting
figure(9)
redChannel = rho' >= 1;greenChannel = rho' < 1;blueChannel=zeros(size(greenChannel));
C = double(cat(3, redChannel, greenChannel, blueChannel));
surf(ks, hns, rho', C)
hold on
w=gobjects(2,1);
w(1)=plot3(nan,nan,nan,'r*');
w(2)=plot3(nan,nan,nan,'g*');
legend(w, ["unstable" "stable"], 'Location','northeast')
xlabel("K_u")
ylabel("hn")
zlabel("\rho")
saveas(gcf,'figures4report/q4-series2.png')




