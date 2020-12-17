% s = CNC_Emulator;


s = serial('/dev/ttyACM0');

fopen(s)

fprintf(s, 'G17 G20 G90 G94 G54')
pause(5)
fprintf(s, 'G1 x0 y0 z.5 F250')
%Raise Cnc
pause(5)
fprintf(s, 'G1 x3.25 y1.75 z.5 F250')
%move cnc above cube
pause(5)
fprintf(s, 'G1 x3.25 y1.75 z4 F250')
%move cnc downward and close arms
pause(5)
fprintf(s, 'G1 x3.25 y1.75 z.5 F250')
%move cnc up
pause(5)
fprintf(s, 'G1 x0 y0 z-8 F250')
%move back to drop off
pause(5)
fprintf(s, 'G1 x0 y0 z4 F250')
%place cube at dropoff zone
pause(5)
fprintf(s, 'G1 x0 y0 z-8 F250')
%REPEAT
pause(5)
fprintf(s, 'G1 x3.25 y-1.75 z.5 F20')
pause(5)
fprintf(s, 'G1 x3.25 y-1.75 z F250')
pause(5)
fprintf(s, 'G1 x3.25 y-1.75 z.5 F250')
pause(5)
fprintf(s, 'G1 x0 y0 z.5 F250')
pause(5)
fprintf(s, 'G1 x0 y0 z3 F250')
pause(5)
fprintf(s, 'G1 x0 y0 z-8 F250')
pause(5)

fclose(s)
