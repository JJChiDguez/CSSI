PrintFile("setup_EC.h", Sprintf("#ifndef SEC\n#define SEC\n") : Overwrite:=true);
l02_a := Floor(var_A / 2);
l01_a := var_A - l02_a;

l02_b := Floor(var_B / 2);
l01_b := var_B - l02_b;

F_p    := GF(var_p);
P_p<t> := PolynomialRing(F_p);

F_q<i> := ext< F_p| t^2 + 1>;
P_q<X> := PolynomialRing(F_q);

a := F_q![0x1, 0x0];
b := F_q![0x0, 0x0];
E_q := EllipticCurve([F_q | a,b]);
load "./magma/balanced_chain_of_d_isogenies.m";

PT := Random(E_q);
while(Order(PT) ne (var_p+1)) do
	PT := Random(E_q);
end while;
PT := (2^(var_A- 1)) * (3^(var_B - 1)) * PT;

KERNEL := [PT]; idx := 1;
while(KERNEL[idx] ne E_q!0) do
	idx +:= 1; KERNEL[idx] := KERNEL[idx - 1] + PT;
end while;

E_q, _ := VelusFormula(KERNEL, E_q);
a := Coefficients(E_q)[4];
b := Coefficients(E_q)[5];

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(a)[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(a)[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(b)[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(b)[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Curve E_0 = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));


// ---------------------------------------------------------------------
P_a := Random(E_q);
while( Order(P_a) ne (var_p+1) ) do
	P_a := Random(E_q);
end while;
P_a := var_f * (3^var_B) * P_a;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(P_a[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(P_a[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(P_a[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(P_a[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point P_a = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// -----------------------
Q_a := Random(E_q);
while true do
	while( (Order(Q_a) ne (var_p+1)) ) do
		Q_a := Random(E_q);
	end while;
	if( (var_f * (3^var_B) * 2^(var_A-1))*Q_a eq 2^(var_A-1)*P_a) then
		Q_a := Random(E_q);
	else
		break;
	end if;
end while;
Q_a := var_f * (3^var_B) * Q_a;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(Q_a[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(Q_a[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(Q_a[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(Q_a[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point Q_a = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// ----------------------- 
b := Random(2);
k := Random(2^(var_A - 1) - 1);

if(b eq 2) then 
	R_a := (2*k)*P_a + Q_a;
else 
	R_a := P_a + (b*2^(var_A-1) + k)*Q_a;
end if;

// ---------------------------------------------------------------------
P_b := Random(E_q);
while( Order(P_b) ne (var_p+1) ) do
	P_b := Random(E_q);
end while;
P_b := var_f * (2^var_A) * P_b;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(P_b[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(P_b[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(P_b[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(P_b[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point P_b = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// -----------------------
Q_b := Random(E_q);
while true do
	while( (Order(Q_b) ne (var_p+1)) ) do
		Q_b := Random(E_q);
	end while;
	if( ((var_f * (2^var_A) * 3^(var_B-1))*Q_b eq -3^(var_B-1)*P_b) or ((var_f * (2^var_A) * 3^(var_B-1))*Q_b eq 3^(var_B-1)*P_b) ) then
		Q_b := Random(E_q);
	else
		break;
	end if;
end while;
Q_b := var_f * (2^var_A) * Q_b;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(Q_b[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(Q_b[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(Q_b[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(Q_b[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point Q_b = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// -----------------------
b := Random(3);
k := Random(3^(var_B - 1) - 1);

if(b eq 2) then 
	R_b := (3*k)*P_b + Q_b;
else 
	R_b := P_b + (b*3^(var_B-1) + k)*Q_b;
end if;
