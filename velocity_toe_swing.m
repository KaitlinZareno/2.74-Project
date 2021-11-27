function drHs = velocity_toe_swing(in1,in2)
%VELOCITY_TOE_SWING
%    DRHS = VELOCITY_TOE_SWING(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    25-Nov-2021 15:02:03

dths = in1(12,:);
dx = in1(7,:);
dy = in1(8,:);
ls = in2(28,:);
ths = in1(6,:);
drHs = [dx+dths.*ls.*cos(ths);dy+dths.*ls.*sin(ths);0.0];
