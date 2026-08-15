// es_eig_chunk.m -- compute a_P(f) for an INTERLEAVED slice of PrimesUpTo(MAXN,F)
// (chunk IDX takes primes IDX, IDX+NCHUNK, IDX+2*NCHUNK, ...), skipping any already in the
// base file, CHECKPOINTING its partial every 50 primes so a timeout loses little.
// A shell driver launches NCHUNK of these; base + partials are merged.
//
// Args: MAXN, NCHUNK, IDX (0-indexed).  Run: magma MAXN:=15000 NCHUNK:=8 IDX:=0 es_eig_chunk.m
SetColumns(0);
AttachSpec("../CHIMP/CHIMP.spec");
AttachSpec("../EichlerShimuraHMF/src/spec");
MAXN := StringToInteger(MAXN); NCHUNK := StringToInteger(NCHUNK); IDX := StringToInteger(IDX);

label := "2.2.8.1-5000.1-j"; DIR := "es_eig/";
f := LMFDBHMF(label);
F := BaseField(Parent(f));
if not assigned f`hecke_eigenvalues then f`hecke_eigenvalues := AssociativeArray(); end if;
// resume: load base file (already-computed primes) so we skip them
base := DIR cat label cat ".txt";
if OpenTest(base, "r") then LoadEigenvalues(~f, base); end if;
have := f`hecke_eigenvalues;      // snapshot of pre-existing keys (don't rewrite these)

primes := PrimesUpTo(MAXN, F);
mine := [ primes[i] : i in [IDX+1..#primes by NCHUNK] | not IsDefined(have, primes[i]) ];
partfile := DIR cat label cat ".part." cat IntegerToString(IDX);

// only-new cache to write to the partial
fnew := f; fnew`hecke_eigenvalues := AssociativeArray();
t0 := Cputime(); c := 0;
for P in mine do
    try fnew`hecke_eigenvalues[P] := HeckeEigenvalue(f, P); catch err ; end try;
    c +:= 1;
    if c mod 50 eq 0 then _ := WriteEigenvalues(fnew, partfile : Overwrite := true); end if;
end for;
n := WriteEigenvalues(fnew, partfile : Overwrite := true);
printf "chunk %o/%o: %o new eigenvalues (of %o mine), %.1o s\n", IDX, NCHUNK, n, #mine, Cputime(t0);
exit;
