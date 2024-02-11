$fn = 64;

electrode_number = 1;
n = electrode_number;

param1 = 10;
L = param1;

param2 = 1;
r = param2;

param3 = 1;
d = param3;

//for (angle = [45:90:315]){
angle = 45 + 90*n;
rotate([0,0,angle])translate([r,0,0])cylinder(h = L, r = d/2,center = true);