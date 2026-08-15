// es_eig_bound.m -- print "BOUND:<n>" where n is the current contiguous norm bound of the
// eigenvalue cache (for PROGRESS.md). Run: magma -b es_eig_bound.m
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");
f := LMFDBHMFwithEigenvalues("2.2.8.1-5000.1-j", "es_eig/");
printf "BOUND:%o\n", NormBoundOnComputedEigenvalues(f);
exit;
