(* Incoming directions *)
a[r_] := Sqrt[1 - r^2];
xi = {
    {1, 1, 0},
    {1, a[s], s},
    {1, a[s], -s}};

(* kappas *)
xi0 = {1, -a[r], r};
ks = LinearSolve[Transpose[xi], xi0];
k1 = ks[[1]];
k2 = ks[[2]];
k3 = ks[[3]];

(* Direct sum decomposition of outgoing direction *)
eta = DiagonalMatrix[{k1, k2, k3}].xi;

(* Incoming 1-form components *)
Y = {
    {0, 0, 1},
    {s, 0, 1},
    {-s, 0, 1}};

(* Signs for Minkowski metric *)
m = {-1, 1, 1};

(* Two-fold interaction coefficients *)
cHalf[k_, l_, b_] :=  Sum[m[[a]]*(
    2 eta[[l, a]] Y[[k, a]] Y[[l, b]] +
    -eta[[l, b]] Y[[k, a]] Y[[l, a]]
    ), {a, 1, 3}];
(* The minus sign corresponds to commutator in the opposite order *)
c[k_, l_, b_] := cHalf[k, l, b] - cHalf[l, k, b];

(* Two-fold interactions *)
p[x_] := -(x[[1]])^2 + (x[[2]])^2 + (x[[3]])^2;
YY12[b_] := c[1, 2, b]/p[eta[[1]] + eta[[2]]];
YY13[b_] := c[1, 3, b]/p[eta[[1]] + eta[[3]]];
YY23[b_] := c[2, 3, b]/p[eta[[2]] + eta[[3]]];

(* Lie algebra components of Y123 \
   Components are indexed by the outer factor x in the nested \
   commutator [x,[y,z]]. \
   The cubic term is discarded as explained in the paper. *)
YYs = {YY23, YY13, YY12};
ees = {eta[[2]] + eta[[3]], eta[[1]] + eta[[3]], eta[[1]] + eta[[2]]};
YYY[x_, b_] := Sum[m[[a]]*(
    -2 eta[[x, a]] Y[[x, b]] YYs[[x]][a] +
    +eta[[x, b]] Y[[x, a]] YYs[[x]][a] +
    +2 ees[[x, a]] Y[[x, a]] YYs[[x]][b] +
    -ees[[x, b]] Y[[x, a]] YYs[[x]][a]
     ), {a, 1, 3}];

(* Print the components *)
pairs = Flatten[Table[{x, b}, {x, 1, 3}, {b, 1, 2}], 1];
printYYY[p_] := 
    Flatten[{p, FullSimplify[Series[Apply[YYY, p], {s, 0, 0}]]}]
TableForm[Map[printYYY, pairs], 
    TableHeadings -> {None, {"x", "beta+1", "Y_{(123)}"}}, 
    TableSpacing -> {5, 2}]

