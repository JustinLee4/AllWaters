load results/test.pdb
hide all
show spheres
set sphere_scale, 0.2
color red, b < .5 and test
color white, b > 0.5 and test
set seq_view, 1, test
load pdbfiles/L-sub_met.pdb
set seq_view, 0, L-sub_met
show seq_view
refresh
save results/test.pse
