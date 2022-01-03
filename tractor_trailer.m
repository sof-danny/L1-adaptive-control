function [dxdt] = tractor_trailer(t,x, V, L1, L2, a, betaf, Vlr, Vsr, Vsi, R1, ...
                                 sigma, zeta, phid, r, omegan, zetad, gamma, gammay, ...
                                 gammar, gammaF, ri, am, bm, taud, P,b, alphax, alphar, alphaf,...
                                 kx, kg, kD, AmL1, bL1, gammaL1, kxi, kgi, kDi, amL1i, gammaiL1, useL1, ...
                                 useprojection)


dxdt = zeros(36,1);

xt     = x(1);
yt     = x(2);
thetat = x(3);
xi     = x(4);
yi     = x(5);
phi    = x(6);

Los     = x(7);
thetaos = x(8);
phios   = x(9);

% reference model for the tractor
xm1 = x(10);
xm2 = x(11);

x1  = x(12); %never actually used
xi2 = x(13); %never actually used

% adaptive tuning laws
L1hat = x(14);

% Trailer MRAC
ym = x(15);
kyhat = x(16);
krhat = x(17);
kFhat = x(18);
y = x(19);

% Tractor Direct MRAC
Mxhat = [x(20), x(21)];
Mrhat = x(22);
Mfhat = x(23);

% Tractor L1 adaptive control
x1tL1hat = x(24); 
x2tL1hat = x(25);
xtL1hat = [x1tL1hat; x2tL1hat];
uad = x(26);
uaddot = x(27);

omegaL1hat = x(28);
thetaL1hat = x(29);
sigmaL1hat = x(30);

% Trailer L1 adaptive control
xiL1hat = x(31);
uadi = x(32);
uadidot = x(33);
omegaiL1hat = x(34);
thetaiL1hat = x(35);
sigmaiL1hat = x(36);


% Tractor control law (includes indirect MRAC, direct MRAC, and L1 adaptive
% control, where the unused laws are commented out
f = (V^2-(sigma(t)*V*sin(thetaos))^2)/(R1(t)+Los); %changed to bring in the disturbances
g = -sigma(t)*V*sqrt(V^2-(sigma(t)*V*sin(thetaos))^2);
e1 = xm1-Los; %Changed to bring in the disturbance
e2 = xm2-(-sigma(t)*V*sin(thetaos)); %changed to bring in the disturbance
mu = 3;
z = e2 + mu*e1;
beta = e1-omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t) + mu*e2;

xtractor = [Los; -sigma(t)*V*sin(thetaos)];
etractor = [e1; e2];

if abs(g) <0.0001
    disp("errors are approaching");
end

xt1L1 = Los; %Measured value of xt1
xt2L1 = -sigma(t)*V*sin(thetaos); %Measured value of xt2
xtL1 = [xt1L1;xt2L1];

%if abs(g) > 0.0001
%u = L1hat/g*(e1-omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t)-f+mu*e2+z+z*mu^2); %(1/1+exp(-10*(abs(g)-1)))*(1/1+exp(-10*(abs(z)-1))) 
%u = Mxhat*xtractor + Mrhat*r(t) + Mfhat;

%Trailer control law
W1 = [-V/L2*cos(phid(t))+sigma(t)*V*a/(R1(t)*L2)*sin(phid(t)); -V/L2*sin(phid(t))-sigma(t)*V*a/(R1(t)*L2)*cos(phid(t)); -sigma(t)*V/R1(t)];
W2 = [V/L2*sin(phid(t))+sigma(t)*V*a/(R1(t)*L2)*cos(phid(t)); -V/L2*cos(phid(t))+sigma(t)*V*a/(R1(t)*L2)*sin(phid(t))];
Phi1 = [sin(thetaos); cos(thetaos); 1];
Phi2 = [sin(thetaos); cos(thetaos)];

if useL1 == 1
    u = uad + kx'*xtL1;
    ui = uadi + kxi*y;
else % Uses Indirect and Direct MRAC
    u = L1hat/g*(e1-omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t)-f+mu*e2+z+z*mu^2);
    ui = kyhat*y+krhat*ri(t)+kFhat;
end

%elseif abs(g) > 0.1
%    u = L1hat/g*(e1-omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t)-f+mu*e2+z);
%else
%u = 0;
%end




% if abs(W2'*Phi2) > 0.01
%     ui = kyhat*y+krhat*ri(t)+kFhat;
% else
%     ui = 0;
% end
%ui = kyhat*y+krhat*ri(t)+kFhat; %Direct MRAC Trailer control law
%ui = uadi + kxi*y; %L1 adaptive control law
 
%Control inputs
%u = u; %+ sin(3*pi*t);
%ui = ui; %+ sin(2*pi*t);
delta = atan(u); %delta  = 0;
deltai = atan(ui);
global u_record t_record adaptive_params_target_record g_rec beta_rec z_rec tractor_params_target_record
u_record = [u_record; delta deltai u ui];
t_record = [t_record; t];
adaptive_params_target_record = [adaptive_params_target_record; (-am)/(W2'*Phi2), (bm)/(W2'*Phi2),(-W1'*Phi1)/(W2'*Phi2)];
g_rec = [g_rec; g];
beta_rec = [beta_rec; beta];
z_rec = [z_rec; z];

tractor_params_target_record = [tractor_params_target_record; -L1/g*omegan^2, -L1/g*2*zeta*omegan, L1/g*omegan^2, -L1/g*f];

%if abs(thetaos)-0.01 > pi/2 || abs(delta)-0.1 > pi/2
%    disp("we approaching a problem")
%end

if t_record(end) > taud && taud > 0
    [M,I] = min(abs(t_record - (t-taud)));
    delta  = u_record(I, 1); %delta  = 0;
    deltai = u_record(I, 2);
end

%Stationary reference frame ODEs
thetai = thetat + phi;
dxdt(1) = (V-Vlr)*cos(thetat)-Vsr*sin(thetat);
dxdt(2) = (V-Vlr)*sin(thetat)+Vsr*cos(thetat);
Thetat  = (V-Vlr)/L1*tan(delta + betaf(t)) + Vsr/L1;
dxdt(3) = Thetat;

N = L2*cos(deltai);
M1 = -(V-Vlr)*sin(deltai+phi);
M2 = Vsr*cos(deltai+phi);
M3 = -Thetat*(a*cos(deltai+phi)+L2*cos(deltai));
M4 = -Vsi;

M3p = -a*Thetat*cos(deltai+phi);

Thetai = 1/N*(M1+M2+M3p+M4);

dxdt(4) = dxdt(1) + Thetat*a*sin(thetat) + Thetai*L2*sin(thetai);
dxdt(5) = dxdt(2) - Thetat*a*cos(thetat) - Thetai*L2*cos(thetai);
dxdt(6) = 1/N*(M1+M2+M3+M4);

% Rotating reference frame ODEs
dxdt(7) = -sigma(t)*abs(V-Vlr)*sin(thetaos)-sigma(t)*zeta*Vsr*cos(thetaos);
dxdt(8) = (V-Vlr)/L1*tan(delta+betaf(t))+Vsr/L1-sigma(t)*abs(V-Vlr)*cos(thetaos)/(R1(t)+Los)+sigma(t)*zeta*Vsr*sin(thetaos)/(R1(t)+Los);
dxdt(9) = 1/(L2*cos(deltai))*( -(V-Vlr)*sin(deltai+phios+phid(t))+Vsr*cos(deltai+phios+phid(t))-((V-Vlr)/L1*tan(delta+betaf(t))+Vsr/L1)*(a*cos(deltai+phios+phid(t))+L2*cos(deltai))-Vsi );



% Vector of tractor time derivatives in the rotating reference frame with
% coordinate change

xt1L1dot = dxdt(7);
xt2L1dot = -sigma(t)*V*dxdt(8)*cos(thetaos);
xtL1dot = [xt1L1dot;xt2L1dot]; %used for higher order filters

% MRAC Tractor 
dxdt(10) = xm2;
dxdt(11) = -omegan^2*xm1-2*zetad*omegan*xm2+omegan^2*r(t);
dxdt(12) = xi2;
dxdt(13) = f+1/L1*g*u; %%%%% This is L1 because what this is simulating is -sigma(t)*V*sin(thetaos), and thetaos is measureable. 


%if L1hat > 0.1 || (L1hat == 0.1 && z*beta+z^2+z^2*mu^2 > 0)
    %dxdt(14) = gamma*Projection(L1hat,z*beta+z^2+z^2*mu^2,6,1); 
dxdt(14) = gamma*z*beta+z^2+z^2*mu^2;
%else
    %dxdt(14) = 0;
%    dxdt(14) = gamma*z*beta+z^2+z^2*mu^2;
%end

% MRAC trailer
dxdt(15) = -am*ym + bm*ri(t);

e = ym - y;

dxdt(16) = gammay*sign(W2'*Phi2)*y*e;
dxdt(17) = gammar*sign(W2'*Phi2)*ri(t)*e;
dxdt(18) = gammaF*sign(W2'*Phi2)*e;

%dxdt(16) = gammay*Projection(kyhat,(-1)*y*e, 3, 0.1);
%dxdt(17) = gammar*Projection(krhat,(-1)*ri(t)*e, 3, 0.1);
%dxdt(18) = gammaF*Projection(kFhat,(-1)*e, 1, 0.1);

dxdt(19) = W1'*Phi1+W2'*Phi2*ui;




dxdt(20:21) = -sigma(t)*(alphax*xtractor*etractor'*P*b)';
dxdt(22) = -sigma(t)*alphar*r(t)*etractor'*P*b;
dxdt(23) = -sigma(t)*alphaf*etractor'*P*b;


%Tractor L1 adaptive controller system of equations
etahat = omegaL1hat*uad + thetaL1hat*max(abs(xtL1)) + sigmaL1hat;
xtL1hatdot = AmL1*xtL1+bL1*(etahat);
dxdt(24) = xtL1hatdot(1);
dxdt(25) = xtL1hatdot(2);

%%%%% HERE is where the D(s) is implemented %%%%%%
rhat = etahat-kg*r(t);
dxdt(26) = -kD*rhat;
dxdt(27) = 0;

eL1 = xtL1hat-xtL1;

if useprojection == 1
    dxdt(28) = Projection(omegaL1hat, -gammaL1*eL1'*P*bL1*uad, 2, 0.1); %omega-hat-tractor
    dxdt(29) = Projection(thetaL1hat,-gammaL1*eL1'*P*bL1*max(abs(xtL1)), 2, 0.1); %theta-hat-tractor
    dxdt(30) = Projection(sigmaL1hat, -gammaL1*eL1'*P*bL1, 2, 0.1); %sigma-hat-tractor
else
    dxdt(28) = -gammaL1*eL1'*P*bL1*uad; %omega-hat-tractor
    dxdt(29) = -gammaL1*eL1'*P*bL1*max(abs(xtL1)); %theta-hat-tractor
    dxdt(30) = -gammaL1*eL1'*P*bL1; %sigma-hat-tractor
end
%Trailer L1 adaptive controller system of equations
etaihat = omegaiL1hat*uadi + thetaiL1hat*abs(y)+sigmaiL1hat;
xiL1hatdot = amL1i*y + etaihat;
dxdt(31) = xiL1hatdot;

%%%%% HERE is where the Di(s) is implemented %%%%%%
rihat = etaihat-kgi*ri(t);
dxdt(32) = -kDi*rihat;
dxdt(33) = 0;

eiL1 = xiL1hat-y;

if useprojection == 1
    dxdt(34) = Projection(omegaiL1hat,-gammaiL1*eiL1*uadi,1,0.1); %omega-hat-trailer
    %dxdt(34) = -sign(W2'*Phi2)*gammaiL1*eiL1*uadi;
    dxdt(35) = Projection(thetaiL1hat,-gammaiL1*eiL1*abs(y),1,0.1); %theta-hat-trailer
    dxdt(36) = Projection(sigmaiL1hat,-gammaiL1*eiL1,1,0.1); %sigma-hat-trailer
else
    dxdt(34) = -gammaiL1*eiL1*uadi; %omega-hat-trailer
    %dxdt(34) = -sign(W2'*Phi2)*gammaiL1*eiL1*uadi;
    dxdt(35) = -gammaiL1*eiL1*abs(y); %theta-hat-trailer
    dxdt(36) = -gammaiL1*eiL1; %sigma-hat-trailer
end



end

