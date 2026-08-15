// goal2_mod5_g_scope.m -- feasibility scope for the descended classical newform g.
// g: weight 2, level N=320000=2^9*5^4, nebentypus chi_8 (quadratic char of Q(sqrt2)),
// degree-4 Hecke field, with A_g = Res_{Q(sqrt2)/Q}(B_f).  We want g's periods, so first
// check whether the modular-symbols space is even constructible here.
// Step 1 (this script): dimensions via formula (fast), then probe ModularSymbols build.
SetColumns(0);
LOGF := "goal2_mod5_g_scope_out.txt";
PrintFile(LOGF, "# g feasibility scope, N=320000, wt 2, chi_8" : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

N := 320000;
// chi_8 = quadratic character of Q(sqrt2): even, conductor 8, chi(3)=chi(5)=-1.
chi8 := KroneckerCharacter(8);
LOG(Sprintf("chi8: modulus %o, order %o, conductor %o, even? %o",
            Modulus(chi8), Order(chi8), Conductor(chi8), IsEven(chi8)));
LOG(Sprintf("  chi8(3)=%o chi8(5)=%o chi8(7)=%o", chi8(3), chi8(5), chi8(7)));

// extend chi8 to modulus N, over the rationals (values +-1 => base field FldRat)
G := DirichletGroup(N, RationalField());
chi := G ! chi8;
LOG(Sprintf("chi mod N: modulus %o, conductor %o, order %o", Modulus(chi), Conductor(chi), Order(chi)));

// dimensions via formula (fast, no space built)
t := Cputime();
dfull := DimensionCuspForms(chi, 2);
dnew  := DimensionNewCuspForms(chi, 2);
LOG(Sprintf("dim S_2(N,chi) FULL cuspidal = %o", dfull));
LOG(Sprintf("dim S_2(N,chi) NEW           = %o   [%.1o s]", dnew, Cputime(t)));

// modular-symbols dimension is ~2x the cuspidal (both +/- signs) unless we fix a sign.
// Probe: try to BUILD ModularSymbols(chi,2, sign 1) and time one Hecke matrix.  This is
// the real feasibility gate (memory/time).  Guarded so the script still logs if it OOMs.
LOG("");
LOG("probing ModularSymbols(chi, 2, +1) construction (may be heavy)...");
t := Cputime();
try
    M := ModularSymbols(chi, 2, 1);
    LOG(Sprintf("  ModularSymbols built: dim %o  [%.1o s]", Dimension(M), Cputime(t)));
    S := CuspidalSubspace(M);
    LOG(Sprintf("  cuspidal subspace dim %o  [%.1o s]", Dimension(S), Cputime(t)));
catch err
    LOG(Sprintf("  BUILD FAILED/too big: %o", err`Object));
end try;
exit;
