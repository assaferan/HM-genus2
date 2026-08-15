// Validate the transcribed dataset: (1) A[l]^ss = 1+chi+sigma structure at good primes
// (catches curve-transcription errors), (2) conductor norm consistency.
SetColumns(0);
load "torsion_data.m";
outf := "validate.out";
PrintFile(outf, "# label  struct(nok/ntot)  condnorm-ok  computed-odd-cond-norm" : Overwrite := true);

for row in torsion_data do
    d := row[1]; l := row[2]; lab := row[3]; fc := row[4]; hc := row[5]; N := row[6]; cc := row[7];
    K<a> := QuadraticField(d); OK := Integers(K); R<x> := PolynomialRing(K); Fl := GF(l);
    f := R![K| c[1]+c[2]*a : c in fc]; hh := R![K| c[1]+c[2]*a : c in hc];
    C := HyperellipticCurve(f, hh);
    // structure check at up to 8 good primes
    nok := 0; ntot := 0;
    for pp in PrimesUpTo(60) do
        if pp eq 2 or pp eq l then continue; end if;
        for tup in Factorization(pp*OK) do
            P := tup[1]; k,red := ResidueClassField(P); Rk := PolynomialRing(k); RT := PolynomialRing(Fl);
            ok := true; Ck := 0;
            try Ck := HyperellipticCurve(Rk![red(c):c in Coefficients(f)], Rk![red(c):c in Coefficients(hh)]);
            catch e; ok := false; end try;
            if not ok then continue; end if;
            wp := RT ! LPolynomial(Ck);
            cp := RT ! [ Coefficient(wp,4-i) : i in [0..4] ];    // char poly of Frob (monic, deg 4)
            Nn := Norm(P) mod l;
            // 1 + chi + sigma: (T-1)(T-Nn) must DIVIDE cp mod l (not merely share the roots --
            // when Nn=1 root-only checks collapse), and the residual quadratic sigma must be the
            // 2-dim piece with det = chi_l, i.e. constant term Nn.
            T := RT.1; q := (T - 1)*(T - Fl!Nn);
            div_ok, sig := IsDivisibleBy(cp, q);
            if div_ok and Degree(sig) eq 2 and Coefficient(sig,0) eq Fl!Nn then nok +:= 1; end if;
            ntot +:= 1;
        end for;
        if ntot ge 8 then break; end if;
    end for;
    // conductor norm
    if d eq 2 then
        Icond := ideal<OK | cc[1]+cc[2]*a>;
        cnorm := Norm(Icond); condok := (cnorm eq N);
    else
        // odd part of the conductor from the curve (Q(sqrt3))
        cnorm := 0; condok := false;
        try
            Cond := Conductor(C);            // full conductor ideal
            odd := 1*OK;
            for pe in Factorization(Cond) do
                if Norm(pe[1]) mod 2 ne 0 then odd := odd * pe[1]^pe[2]; end if;
            end for;
            cnorm := Norm(odd); condok := (cnorm eq N);
        catch e; cnorm := -1; end try;
    end if;
    line := Sprintf("%-11o  struct %o/%o  cond-ok=%o  cnorm=%o (stated %o)", lab, nok, ntot, condok, cnorm, N);
    printf "%o\n", line; PrintFile(outf, line);
end for;
PrintFile(outf, "# done"); printf "DONE\n"; quit;
