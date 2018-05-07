PrintFile("setup_EC.h", "//---------------------------------------------");
// ---------------------------------------------------------------------
E_B := d_pow_e_isogenous_curve(R_b, E_q, 3, var_B);

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(Coefficients(E_B)[4])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(Coefficients(E_B)[4])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(Coefficients(E_B)[5])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(Coefficients(E_B)[5])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Curve E_B = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

psiP_a, psiQ_a := double_eval_d_pow_e_isogeny(R_b, P_a, Q_a, E_q, 3, var_B);

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(psiP_a[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(psiP_a[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(psiP_a[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(psiP_a[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point psi_P_a = { { {%h, %h}, {%h, %h} },\n\t\t\t { {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(psiQ_a[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(psiQ_a[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(psiQ_a[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(psiQ_a[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point psi_Q_a = { { {%h, %h}, {%h, %h} },\n\t\t\t { {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// -----------------------
S_b := Random(E_B);
while( Order(S_b) ne (var_p+1) ) do
	S_b := Random(E_B);
end while;
S_b := var_f * (2^var_A) * S_b;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(S_b[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(S_b[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(S_b[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(S_b[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point S_b = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));

// -----------------------
T_b := Random(E_B);
while true do
	while( (Order(T_b) ne (var_p+1)) ) do
		T_b := Random(E_B);
	end while;
	if( ((var_f * (2^var_A) * 3^(var_B-1))*T_b eq -3^(var_B-1)*S_b) or ((var_f * (2^var_A) * 3^(var_B-1))*T_b eq 3^(var_B-1)*S_b) ) then
		T_b := Random(E_B);
	else
		break;
	end if;
end while;
T_b := var_f * (2^var_A) * T_b;

//---
HEX_AX := Intseq(IntegerRing()!Eltseq(T_b[1])[1], 2^64);
HEX_AX cat:= [0 : i in [1 .. (var_words - #HEX_AX)]];
HEX_AY := Intseq(IntegerRing()!Eltseq(T_b[1])[2], 2^64);
HEX_AY cat:= [0 : i in [1 .. (var_words - #HEX_AY)]];

HEX_BX := Intseq(IntegerRing()!Eltseq(T_b[2])[1], 2^64);
HEX_BX cat:= [0 : i in [1 .. (var_words - #HEX_BX)]];
HEX_BY := Intseq(IntegerRing()!Eltseq(T_b[2])[2], 2^64);
HEX_BY cat:= [0 : i in [1 .. (var_words - #HEX_BY)]];
//---

PrintFile("setup_EC.h", Sprintf("static Point T_b = { \t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {%h, %h}, {%h, %h} },\n\t\t\t{ {0x1, 0x0}, {0x0, 0x0} }\n\t\t   };",
HEX_AX[1], HEX_AX[2], HEX_AY[1], HEX_AY[2], HEX_BX[1], HEX_BX[2], HEX_BY[1], HEX_BY[2]));
PrintFile("setup_EC.h", Sprintf("#endif"));
