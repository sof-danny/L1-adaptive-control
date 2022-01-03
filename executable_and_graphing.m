% For your reference

% t is an array of timestamps for the simulation outputs
% y(:,1) is the x position of the tractor
% y(:,2) is the y position of the tractor
% y(:,3) is the angle of the tractor (angle along the body)
% y(:,4) is the x position of the trailer
% y(:,5) is the y position of the trailer
% y(:,6) is the angle the trailer hitch makes with the tractor
% 
% y(:,7) is the distance between the tractor's center (xt,yt) and the circle it is trying to track (L-offset)
% y(:,8) is the difference between the tractor's angle and the ideal angle for the circle it is trying to track (theta offset)
% y(:,9) is the difference between the trailer hitches' angle and the ideal angle for the circle it is trying to track (phi offset)
% 
% y(:,10) is the model reference for y(:,7) which is tracking r(t)
% y(:,11) is the model reference for y(:,8) which is going to 0
% 
% y(:,12) is the same as y(:,7) if things are going correctly (in the code, it is calculated a slightly different way, so agreement is a debugging tool)
% y(:,13) is -sigma*V*sin(y(:,11)) if things are going correctly (again, a debugging tool of sorts)
% 
% y(:,14) is L1hat, i.e. the estimate of L1 achieved through adaptation
% 
% y(:,15) is the model reference for y(:,9) which is going to ri(t)
% 
% y(:,16) is kyhat i.e. the estimate of the ideal gain for one of the terms in the trailer's controller (The physical meaning of ky is confusing and not-intuitive)
% y(:,17) is krhat i.e. the estimate of the ideal gain for another of the terms in the trailer's controller (The physical meaning of kr is confusing)
% y(:,18) is kFhat (same as above, this is a weird value which appears in the trailer's controller)
% -- All of these terms are estimates of something which when accurate results in cancellation of nonlinearities and decent control
% -- The estimates adapt when ri(t) is a persistent excitation
% 
% y(:,19) is the same as y(:,9) if we assume that the tractor is perfectly following its intended trajectory
% 
% y(:,20:23) are the adaptive parameters for direct MRAC of the tractor.
% They are assuming that the disturbances are zero. I'm not sure what
% happens if we are wrong about this. 

% ALSO I used some global variables to extract the controller values and the values which kyhat, krhat, and kFhat are approximating (They have formulas, those formulas are ugly)
% These are stored in the global variables u_record, t_record, and adaptive_params_target_record. Because the "t" values spit out by ode45() are different than those used by the simulation, 
% t and t_record are not the same!!!!

% Model parameters
V     = 1;
L1    = 3.5;
L2    = 1;
a     = 0.5;

% Disturbance parameters
betaf =@(t) 0;%0.2*sin(t); %*exp(-0.01*t);
Vlr   = 0;
Vsr   = 0;
Vsi   = 0;


% Trajectory path following
R1 = @(t) 15; %15 + 1000*0.5.*(square(t*(2*pi)/40)+1)
sigma = @(t) -1; %1.*(t<35) -1.*(t>=35); % +1 for CCW, -1 for CW

r  = @(t) 1+1.*(t>25); %0.01*randi([0 1]); %(0.1*sin(t) + 0.1*sin(2*t) + 0.1*sin(3*t)).*(t<50); %1.5.*(t<50)+0.5.*(t<70); %1+sin(t); %sin(t); %(0.1*sin(t) + 0.1*sin(2*t) + 0.1*sin(3*t)).*(t<50) ;
rdot = @(t) 0;
ri = @(t) 0.*(t>0); %0.01*randi([0 1]); %(0.1*sin(t) + 0.1*sin(2*t) + 0.1*sin(3*t)).*(t<3000);

zeta = 1;
phid =@(t) -sigma(t)*(pi - atan(R1(t)/a)-acos( (L2^2+a^2)/(2*R1(t)*sqrt(a^2+R1(t)^2)) ));

% Initial Conditions
xt0 = 0;
yt0 = 0;
thetat0 = 0;
phi0 = 0;
xi0 = xt0 - a*cos(thetat0) - L2*cos(thetat0+phi0);
yi0 = yt0 - a*sin(thetat0) - L2*sin(thetat0+phi0);

Los0     = 0;
thetaos0 = 0;
phios0   = phi0-phid(0);

L1hat0 = 2;

ym0 = phios0;
kyhat0 = 0;
krhat0 = 0;
kFhat0 = 0;
y0 = phios0;

Mx1hat0 = -1;
Mx2hat0 = -1.4;
Mrhat0  = 1;
Mfhat0  = -1;

% L1 adaptive state equations
x1tL1hat0 = Los0; 
x2tL1hat0 = -sigma(0)*V*sin(thetaos0);
uad0 = 0;
uaddot0 = 0;

omegaL1hat0 = 0;
thetaL1hat0 = 0;
sigmaL1hat0 = 0;

% Trailer L1 adaptive state equations
xiL1hat0 = phios0; 
uadi0 = 0; 
uadidot0 = 0; 
omegaiL1hat0 = 0; 
thetaiL1hat0 = 0; 
sigmaiL1hat0 = 0;



%Initial conditions
ic = [xt0; yt0; thetat0; xi0; yi0; phi0; Los0; thetaos0; phios0; Los0; ...
      thetaos0; Los0; -sigma(0)*V*sin(thetaos0); L1hat0; ym0; kyhat0; krhat0; ...
      kFhat0; y0; Mx1hat0; Mx2hat0; Mrhat0; Mfhat0; ...
      x1tL1hat0; x2tL1hat0; uad0; uaddot0; omegaL1hat0; thetaL1hat0; ...
      sigmaL1hat0; ...
      xiL1hat0; uadi0; uadidot0; omegaiL1hat0; thetaiL1hat0; sigmaiL1hat0];

% control laws

omegan = 1;
zetad = 0.7;

% Tuning law
gamma = 100;

%Trailer control parameters
gammay = 50;
gammar = 50;
gammaF = 50;

%Trailer MRAC
am = 1;
bm = 1;

% Tractor MRAC
Am = [0 1; -omegan^2 -2*zetad*omegan];
Q = eye(2);
P = lyap(Am',Q);
b = [0;1];

alphax = 20;
alphar = 20;
alphaf = 20;

% L1 adaptive parameters

kx = [omegan^2; 2*zetad*omegan];
cL1 = [1 0];
AmL1 = [0 1; -omegan^2, -2*zetad*omegan];
bL1 = [0;1];
kg = -(cL1*(AmL1\bL1))^-1;

kD = 10; % Tractor controller cut off frequency (>10)
gammaL1 = 10000; %limited to around 10^4, otherwise high frequency noise leaks into the control law. Maybe projection helps with this?

% L1 trailer adaptive parameters
amL1i = -1;
kxi = -amL1i;
kgi = -amL1i;
kDi = 0.1; % Trailer control signal cut off frequency (<10)
gammaiL1 = 1000;


%% Simulation
global u_record t_record adaptive_params_target_record g_rec beta_rec z_rec tractor_params_target_record
u_record = [0, 0, 0, 0];
t_record = [0];
adaptive_params_target_record = [0,0,0];
tractor_params_target_record = [0,0,0,0];
g_rec = [0];
beta_rec = [0];
z_rec = [0];

tspan = [0,20]

taud = 0.15; %0.105 maximum for high gain 100, 50, 50, 50, 0.300 maximum for low gain
useL1 = 1;
useprojection = 1;

opts = odeset('MaxStep',1e-2);
[t,y] = ode45(@(t,x) tractor_trailer(t,x, V, L1, L2, a, betaf, Vlr, ...
                                     Vsr, Vsi, R1, sigma, zeta, phid, r, ...
                                     omegan, zetad, gamma, gammay, gammar, ...
                                     gammaF, ri, am, bm, taud,P,b, alphax, ...
                                     alphar, alphaf, kx, kg, kD, AmL1, bL1, gammaL1,...
                                     kxi, kgi, kDi, amL1i, gammaiL1, useL1, useprojection), tspan, ic, opts);

disp("Simulation completed successfully");

Los     = y(:,7);
thetaos = y(:,8);
phios   = y(:,9);

%% Optional code for reducing the size of the output matrices
%y = y((1:10:end), :);
%t = t(1:10:end);

%% Graphing  stationary reference frame
close all
figure(1)
sn = 10; sn2 = 200;sn3 = 30000;
figure(1)
plot(y(:,1),y(:,2),y(:,4),y(:,5))
% xlim([-10,10]);
% ylim([-10,10]);
hold on
set(gca, 'FontSize', 14)

%plot(R1(0)*cos(0:0.1:2*pi),R1(0)*sin(0:0.1:2*pi)-R1(0),"r--")
plot(R1(0)*cos(0:0.1:2*pi),R1(0)*sin(0:0.1:2*pi)-R1(0),"r--") 

pt  = [y(:,1), y(:,2)];
ptf = [y(:,1)+L1*cos(y(:,3)), y(:,2)+L1*sin(y(:,3))];
ptr = [y(:,1)-a*cos(y(:,3)), y(:,2)-a*sin(y(:,3))];
pti  = [y(:,4),y(:,5)];

px = [ptf(:,1), pt(:,1), ptr(:,1), pti(:,1)];
py = [ptf(:,2), pt(:,2), ptr(:,2), pti(:,2)];

pxt = [ptf(:,1), pt(:,1), ptr(:,1)];
pyt = [ptf(:,2), pt(:,2), ptr(:,2)];

pxi = [ptr(:,1), pti(:,1)];
pyi = [ptr(:,2), pti(:,2)];


if length(t) > sn3
plot(pxt(sn,:),pyt(sn,:),'k.-','LineWidth',2,'MarkerSize',20)
plot(pxi(sn,:),pyi(sn,:),'b.-','LineWidth',2,'MarkerSize',20)

plot(pxt(sn2,:),pyt(sn2,:),'k.-','LineWidth',2,'MarkerSize',20)
plot(pxi(sn2,:),pyi(sn2,:),'b.-','LineWidth',2,'MarkerSize',20)

plot(pxt(sn3,:),pyt(sn3,:),'k.-','LineWidth',2,'MarkerSize',20)
plot(pxi(sn3,:),pyi(sn3,:),'b.-','LineWidth',2,'MarkerSize',20)
end

legend("Tractor", "Trailer", "Path")

xlabel("x (m)")
ylabel("y (m)")
grid on

%% Offset Model graphs
figure(2)
subplot(3,1,1)
plot(t, y(:,7))
hold on
fplot(r, tspan)
grid on
%xlabel("Time (s)")
%ylabel({"Tractor offset","distance (m)"})
%legend("L_{os}","r(t)")
%set(gca, 'FontSize', 14)

subplot(3,1,2)
plot(t, y(:,8))
hold on
fplot(@(t) 0, tspan);
grid on
%xlabel("Time (s)")
%ylabel({"Tractor offset","angle (rad)"})
%legend("\theta_{os}", "y = 0")
%set(gca, 'FontSize', 14)

subplot(3,1,3)
plot(t, y(:,9))
hold on
fplot(ri, tspan)
grid on
%xlabel("Time (s)")
%ylabel({"Trailer offset","angle (rad)"})
%legend("\phi_{os}", "r_{i}(t)")
%set(gca, 'FontSize', 14)

%% Trailer Adaptation
figure(5)
% subplot(3,1,1)
% plot(t, y(:,14))
% hold on; grid on;
% plot(t, t*0+L1);
% xlabel("Time (s)")
% ylabel("Estimated Tractor Length (m)")
% legend("$$\hat{L_{1}}$$","$$L_{1}$$", "interpreter", "LaTeX")
% set(gca, 'FontSize', 14)


%subplot(3,1,[2,3])
plot(t_record,adaptive_params_target_record(:,1), "r--"); hold on; grid on;
plot(t_record,adaptive_params_target_record(:,2), "b--")
plot(t_record,adaptive_params_target_record(:,3), "g--")
plot(t, y(:,16),"r");
plot(t, y(:,17),"b");
plot(t, y(:,18),"g");
xlabel("Time (s)")
ylabel("Parameter estimates (-/-)")


legend("$$k_{y}^*$$","$$k_{r}^*$$","$$k_{F}^*$$", "$$\hat{ky}$$", "$$\hat{kr}$$", "$$\hat{kF}$$", "interpreter", "LaTeX")
set(gca, 'FontSize', 14)
% 
% %%
% figure(6)
% plot(t, y(:,14))
% hold on; grid on;
% plot(t, t*0+L1);
% xlabel("Time (s)")
% ylabel("Estimated Tractor Length (m)")
% legend("$$\hat{L_{1}}$$","$$L_{1}$$", "interpreter", "LaTeX")
% set(gca, 'FontSize', 14)

%% Tractor Adaptation
Mxhatstar = tractor_params_target_record(:,1:2);
Mrhatstar = tractor_params_target_record(:,3);
Mfhatstar = tractor_params_target_record(:,4);


figure(6)

subplot(3,1,1)
plot(t_record,Mxhatstar(:,1), "r--"); hold on; grid on;
plot(t_record,Mxhatstar(:,2), "b--");
plot(t, y(:,20),"r");
plot(t, y(:,21),"b");
xlabel("Time (s)")
ylabel("Parameter estimates (-/-)")
ylim([-10,10])


subplot(3,1,2)
plot(t_record,Mrhatstar, "g--"); hold on; grid on;
plot(t, y(:,22),"g");
xlabel("Time (s)")
ylabel("Parameter estimates (-/-)")
ylim([-10,10])

subplot(3,1,3)
plot(t_record,Mfhatstar, "k--"); hold on; grid on;
plot(t, y(:,23),"k");
xlabel("Time (s)")
ylabel("Parameter estimates (-/-)")
ylim([-10,10])


%set(gca, 'FontSize', 14)

%% Graphing trailer and tractor adaptive L1 parameters
figure(88)
subplot(2,1,1)
plot(t, y(:, 28))
hold on
plot(t, y(:,29))
plot(t, y(:,30))
legend("omegat","thetat", "sigmat")

subplot(2,1,2)
plot(t, y(:, 34))
hold on
plot(t, y(:,35))
plot(t, y(:,36))
legend("omegai","thetai", "sigmai")

%% Legacy graphing tools
% figure(7)
% subplot(4,1,1)
% plot(phios_05{1},phios_05{2}); hold on;xlim([0,20]); ylim([-0.5, 1.5])
% xlabel("Time (s)")
% ylabel("\phi_{os}")
% legend("r_{i}(t) = 0.5")
% subplot(4,1,2)
% plot(phios_10{1},phios_10{2});xlim([0,20]); ylim([-0.5, 1.5])
% xlabel("Time (s)")
% ylabel("\phi_{os}")
% legend("r_{i}(t) = 1.0")
% subplot(4,1,3)
% plot(phios_15{1},phios_15{2});xlim([0,20]); ylim([-0.5, 1.5])
% xlabel("Time (s)")
% ylabel("\phi_{os}")
% legend("r_{i}(t) = 1.5")
% subplot(4,1,4)
% plot(phios_20{1},phios_20{2});xlim([0,20]); ylim([-0.5, 1.5])
% xlabel("Time (s)")
% ylabel("\phi_{os}")
% legend("r_{i}(t) = 2.0")



% %%
% figure(3)
% subplot(3,1,1)
% plot(t, y(:,14))
% hold on; grid on;
% plot(t, t*0+L1);
% xlabel("Time (s)")
% ylabel("Estimated Tractor Length (m)")
% legend("$$\hat{L_{1}}$$","$$L_{1}$$", "interpreter", "LaTeX")
% set(gca, 'FontSize', 14)
% subplot(3,1,[2,3])
% plot(t_record(1:1000:end),adaptive_params_target_record(1:1000:end,1), "r--"); hold on; grid on;
% plot(t_record(1:1000:end),adaptive_params_target_record(1:1000:end,2), "b--")
% plot(t_record(1:1000:end),adaptive_params_target_record(1:1000:end,3), "g--")
% plot(t, y(:,16),"r");
% plot(t, y(:,17),"b");
% plot(t, y(:,18),"g");
% xlabel("Time (s)")
% ylabel("Parameter estimates (-/-)")
% 
% legend("$$k_{y}^*$$","$$k_{r}^*$$","$$k_{F}^*$$", "$$\hat{ky}$$", "$$\hat{kr}$$", "$$\hat{kF}$$", "interpreter", "LaTeX")
% set(gca, 'FontSize', 14)
% %%
% figure(4)
% plot(t_record, rad2deg(u_record(:,1)))
% hold on; grid on
% plot(t_record, rad2deg(u_record(:,2)))
% legend("\delta","\delta_{i}","location","best")
% ylabel("Wheel angle ("+char(176)+")")
% xlabel("Time (s)")
% set(gca, 'FontSize', 14)
% % %% figure(2)
% % 
% % subplot(7,2,2);
% % plot(t, Los); hold on;
% % plot(t, r(t));
% % legend("Los", "r(t)")
% % subplot(7,2,4);
% % plot(t, thetaos, t, phios)
% % hold on 
% % plot(t, r(t), "--")
% % legend("Thetaos", "\phi_{os}")
% % 
% % subplot(7,2,6)
% % plot(t,y(:,14),t, ones(size(t))*L1);
% % legend("L1hat", "L1");
% % 
% % subplot(7,2,7:10)
% % plot(t_record,adaptive_params_target_record(:,1), "r--"); hold on; grid on;
% % plot(t_record,adaptive_params_target_record(:,2), "b--")
% % plot(t_record,adaptive_params_target_record(:,3), "g--")
% % plot(t, y(:,16),"r");
% % plot(t, y(:,17),"b");
% % plot(t, y(:,18),"g");
% % 
% % 
% % legend("ky^*","kr^*","kF^*", "kyhat", "krhat", "kFhat")
% % 
% % subplot(7,2,11)
% % plot(t_record, u_record(:,1));
% % hold on;
% % plot(t_record, u_record(:,2));
% % legend("Front wheel angle", "Trailer Wheel angle")
% % 
% % 
% % subplot(7,2,12)
% % plot(t, y(:,15)); 
% % hold on; 
% % plot(t, y(:,19));
% % legend("ym", "y");
% % 
% % % figure(2)
% % % subplot(3,1,1)
% % % plot(t, -sigma(t).*V.*sin(y(:,8)),t, y(:,13))
% % % subplot(3,1,2)
% % % plot(t, y(:,7),t, y(:,12))
% % % subplot(3,1,3)
% % % plot(t,y(:,14))
% % 
% % % Graphing interesting values
% % 
% % %%
% % %figure(2)
% % %plot(t, R1(t))
% 
% 