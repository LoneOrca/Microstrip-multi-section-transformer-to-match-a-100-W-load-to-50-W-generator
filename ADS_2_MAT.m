file = 'Micro_V1.s2p';
file2  = 'Microstep_MSTEP_V2.s2p';
g = read(rfdata.data,file);
h = read(rfdata.data,file2);
figure
plot(g,'s11','s21','db');
hold all
plot(h,'s11','s21','db');