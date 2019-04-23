%This script models the dynamics of the Omni-directional mobile robot
%with three wheeled traingular arrangment 

%% Parameters definitions
m = 1.5457;    %robot mass [kg]
Iz = 13493290.02;   %robot moment of inertia [kg m^2]
L = 1;          %radius of robot body [m]
R = 1;          %wheel radius
n = 64;         %gear ratio
Ra = 1;         %armature resistance
k2 = 1;         %motor torque constant
k3 = 1.4;       %back emf constant
J0 = 5.7E-7;    %combined inertia of the motor, gear train and wheel referred to the motor shaft
b0 = 1;         %viscous-friction coefficient of the motor, gear and wheel combination

H = [1/m, 0, 0;
    0, 1/m, 0;
    0, 0, 1/Iz]; %mass matrix

B = [0, cos(pi/6), -cos(pi/6);...
    -1, sin(pi/6), sin(pi/6);...
    L, L, L];      %this matrix relates the 3 forces on the wheels to the body-fixed accelerations

G = eye(3)+H*B*B'*(n^2*J0/R^2);           %just a lumped matrix for calculations
Ginv = inv(G);
A1 = Ginv*H*B*B'*(k2*k3/Ra+b0)*(n^2/R^2); %lumped matrix for easy computations
A2 = Ginv*H*B*(k2*n/R/Ra);                %lumped matrix for easy computations

%% Numerical simulation of the EOM

%states definition
%u -> x(1)
%v -> x(2)
%r -> x(3)
%x -> x(4)
%y -> x(5)
%psi -> x(6)


tsim = 10;            %simulation time (s)
timespan = [0 tsim];
IC = [1;1;1;0;0;0];   %initial conditions

%control inputs
tcontrol = 0:.01:tsim;
for i = 1:length(tcontrol)
    E1(i) = 0;   %applied voltage on motor1
    E2(i) = 0;   %applied voltage on motor2
    E3(i) = 0;   %applied voltage on motor3
end

gvEOM = @(t,x) [Ginv(1,1)*x(3)*x(2)-Ginv(1,2)*x(3)*x(1)-A1(1,1)*x(1)-A1(1,2)*x(2)-A1(1,3)*x(3)+...
    A2(1,1)*interp1(tcontrol,E1,t)+A2(1,2)*interp1(tcontrol,E2,t)+A2(1,3)*interp1(tcontrol,E3,t);...
    Ginv(2,1)*x(3)*x(2)-Ginv(2,2)*x(3)*x(1)-A1(2,1)*x(1)-A1(2,2)*x(2)-A1(2,3)*x(3)+...
    A2(2,1)*interp1(tcontrol,E1,t)+A2(2,2)*interp1(tcontrol,E2,t)+A2(2,3)*interp1(tcontrol,E3,t);...
    Ginv(3,1)*x(3)*x(2)-Ginv(3,2)*x(3)*x(1)-A1(3,1)*x(1)-A1(3,2)*x(2)-A1(3,3)*x(3)+...
    A2(3,1)*interp1(tcontrol,E1,t)+A2(3,2)*interp1(tcontrol,E2,t)+A2(3,3)*interp1(tcontrol,E3,t);...
    cos(x(6))*x(1)-sin(x(6))*x(2);...
    sin(x(6))*x(1)+cos(x(6))*x(2);...
    x(3)];

[T ,Y] = ode45( gvEOM,timespan ,IC);

%% Plotting the results
plot(T,Y(:,1),T,Y(:,2),T,Y(:,3),T,Y(:,4),T,Y(:,5),T,Y(:,6))
legend('u','v','r','x','y','psi')
