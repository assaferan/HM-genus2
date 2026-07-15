// goal2_FA_recognize.m -- high-precision attempt to recognize the degree-400 F_A
// resolvent R(X) in Q[X] (idx 33), with RESUME checkpointing (survives host reboots).
//
// Invariant matches goal2_FA_spike.m: theta_line = sum over the 3 Kummer points
// {[P],[2P],[3P]} of t = s1 + rho*s2  (rho=3), s1,s2 = Mumford u-poly of the divisor.
// Each theta_line is checkpointed (integer mantissa + base-10 exponent, exact) as it
// is computed; a restart reloads the checkpoint and resumes from the first missing line.
//
// Run: magma prc:=1200 goal2_FA_recognize.m
SetColumns(0);
prec := StringToInteger(prc);
rho  := 3;
SD   := prec + 10;                                   // stored significant digits
QQ := Rationals(); Px<x> := PolynomialRing(QQ);
C := HyperellipticCurve(-10*x^5 + 15*x^4 - 30*x^3 + 10*x^2 + 60*x + 9);
LOGF := "goal2_FA_recognize_out.txt";
CKPT := Sprintf("goal2_FA_ckpt_p%o.txt", prec);
PrintFile(LOGF, Sprintf("# F_A degree-400 resolvent recognition, prec=%o", prec) : Overwrite := true);
procedure LOG(s) printf "%o\n", s; PrintFile(LOGF, s); end procedure;

t0 := Cputime();
AJ  := AnalyticJacobian(C : Precision := prec);
BPM := BigPeriodMatrix(AJ); CC := BaseRing(BPM);
lat := [ Matrix(CC,2,1,[BPM[1,j],BPM[2,j]]) : j in [1..4] ];
LOG(Sprintf("period matrix built (%.1o s)", Cputime(t0)));

// exact integer (mantissa,exponent) codec for full-precision reals
function EncReal(v)
    if v eq 0 then return "0 0"; end if;
    E := Floor(Log(10, Abs(v)));
    M := Round( (v / (CC!10)^E) * (CC!10)^SD );
    return Sprintf("%o %o", M, E);
end function;
function DecReal(mm, ee) return (CC!mm) * (CC!10)^(ee - SD); end function;

function KummerT(w)
    z := (1/7)*&+[ (CC!w[j])*lat[j] : j in [1..4] ];
    pts := FromAnalyticJacobian(z, AJ);
    error if #pts ne 2, "anomalous divisor (#pts <> 2) at w =", w;
    x1 := pts[1][1]; x2 := pts[2][1];
    return (x1+x2) + rho*(x1*x2);                    // s1 + rho*s2
end function;

// enumerate 400 lines (leading nonzero entry = 1)
lines := [];
for v1,v2,v3,v4 in [0..6] do
    w := [v1,v2,v3,v4]; lead := 0;
    for j in [1..4] do if w[j] ne 0 then lead := w[j]; break; end if; end for;
    if lead eq 1 then Append(~lines, w); end if;
end for;

// ---- load checkpoint ----
done := AssociativeArray();          // li -> theta (CC)
try
    fh := Open(CKPT, "r");
    while true do
        s := Gets(fh);
        if IsEof(s) then break; end if;
        p := Split(s, "|");          // "li|Mre Ere|Mim Eim"
        li := StringToInteger(p[1]);
        pr := Split(p[2], " "); pi := Split(p[3], " ");
        done[li] := DecReal(StringToInteger(pr[1]),StringToInteger(pr[2]))
                  + CC.1*DecReal(StringToInteger(pi[1]),StringToInteger(pi[2]));
    end while;
    delete fh;
catch e; end try;
LOG(Sprintf("resume: %o/%o lines already in checkpoint", #Keys(done), #lines));

// ---- compute theta_line for the remaining lines, checkpointing each ----
tt := Cputime();
for li in [1..#lines] do
    if IsDefined(done, li) then continue; end if;
    w := lines[li];
    t := KummerT([w[j] mod 7 : j in [1..4]])
       + KummerT([(2*w[j]) mod 7 : j in [1..4]])
       + KummerT([(3*w[j]) mod 7 : j in [1..4]]);
    done[li] := t;
    PrintFile(CKPT, Sprintf("%o|%o|%o", li, EncReal(Real(t)), EncReal(Imaginary(t))));
    if li mod 25 eq 0 then LOG(Sprintf("  ... line %o/%o (%.0o s, %o done)", li, #lines, Cputime(tt), #Keys(done))); end if;
end for;
thetas := [ done[li] : li in [1..#lines] ];
LOG(Sprintf("all %o theta_line ready", #thetas));

// build resolvent numerically
CCx<X> := PolynomialRing(CC);
R := &*[ X - t : t in thetas ];
LOG(Sprintf("resolvent degree %o built", Degree(R)));

// recognition: smallest height bound at which each coefficient's best rational locks on
bmax := (prec div 2) - 50;
lockthresh := 10^(-(Ceiling(prec*8/10)));   // rational, NOT complex (real<complex is unordered)
recognized := 0; heights := []; failed := []; tracelocked := "no (height > 10^" cat IntegerToString(bmax) cat ")";
for i in [0..Degree(R)] do
    re := Real(Coefficient(R, i));
    locked := false;
    for b in [50*j : j in [1..bmax div 50]] do
        q := BestApproximation(re, 10^b);
        if Abs(re - q) lt lockthresh then
            h := Ceiling(Log(10, Max(Abs(Numerator(q)), Abs(Denominator(q)))+1));
            Append(~heights, h); recognized +:= 1; locked := true;
            if i eq Degree(R)-1 then tracelocked := Sprintf("height 10^%o", h); end if;
            break;
        end if;
    end for;
    if not locked then Append(~failed, i); end if;
    if i mod 40 eq 0 then LOG(Sprintf("  recog %o/%o (%o locked)", i, Degree(R)+1, recognized)); end if;
end for;

LOG("\n=== RECOGNITION RESULT ===");
LOG(Sprintf("coefficients recognized: %o / %o", recognized, Degree(R)+1));
if #heights gt 0 then LOG(Sprintf("recognized heights: min 10^%o, max 10^%o", Min(heights), Max(heights))); end if;
LOG(Sprintf("trace (e_1): %o", tracelocked));
LOG(Sprintf("unrecognized (height > 10^%o): %o coefficients", bmax, #failed));
if recognized eq Degree(R)+1 then
    Rq := Px![ BestApproximation(Real(Coefficient(R,i)), 10^bmax) : i in [0..Degree(R)] ];
    PrintFile("goal2_FA_resolvent.txt", Sprintf("%o", Rq) : Overwrite := true);
    LOG("FULLY RECOGNIZED -- R(X) written to goal2_FA_resolvent.txt");
else
    LOG("NOT fully recognized at this precision (expected -- quantifies the height wall).");
end if;
LOG("DONE.");
exit;
