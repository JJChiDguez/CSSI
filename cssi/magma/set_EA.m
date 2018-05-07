PrintFile("setup_EC.h", "//---------------------------------------------");
// ---------------------------------------------------------------------
E_A := d_pow_e_isogenous_curve(R_a, E_q, 2, var_A);

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(Coefficients(E_A)[4])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(Coefficients(E_A)[4])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(Coefficients(E_A)[5])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(Coefficients(E_A)[5])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Curve E_A = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

phiP_b, phiQ_b := double_eval_d_pow_e_isogeny(R_a, P_b, Q_b, E_q, 2, var_A);

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(phiP_b[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(phiP_b[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(phiP_b[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(phiP_b[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point phi_P_b = { { {%h, %h}, {%h, %h} },\n\t\t\t { {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(phiQ_b[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(phiQ_b[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(phiQ_b[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(phiQ_b[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point phi_Q_b = { { {%h, %h}, {%h, %h} },\n\t\t\t { {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// -----------------------
S_a := Random(E_A);
while( Order(S_a) ne (var_p+1) ) do
	S_a := Random(E_A);
end while;
S_a := var_f * (3^var_B) * S_a;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(S_a[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(S_a[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(S_a[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(S_a[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point S_a = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// -----------------------
T_a := Random(E_A);
while true do
	while( (Order(T_a) ne (var_p+1)) ) do
		T_a := Random(E_A);
	end while;
	if( (var_f * (3^var_B) * 2^(var_A-1))*T_a eq 2^(var_A-1)*S_a) then
		T_a := Random(E_A);
	else
		break;
	end if;
end while;
T_a := var_f * (3^var_B) * T_a;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(T_a[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(T_a[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(T_a[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(T_a[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point T_a = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));
