PrintFile("setup_FF.h", Sprintf("#ifndef SP\n#define SP\n") : Overwrite:=true);
PrintFile("setup_FF.h", Sprintf("static uint8_t expn_a = %o;", var_A));
PrintFile("setup_FF.h", Sprintf("static uint8_t expn_b = %o;\n", var_B));
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
var_p := 2^var_A * 3^var_B * var_f - 1;
//print Intseq(var_p, 2^64):Hex;
PrintFile("setup_FF.h", Sprintf("static uint64_t\t\t  p[] = {%h, %h};", Intseq(var_p, 2^64)[1], Intseq(var_p, 2^64)[2]));
var_r := (2^64)^var_words;
var_R := var_r mod var_p;
neg_R := (-var_R) mod var_p;
//print Intseq(var_R, 2^64):Hex;
PrintFile("setup_FF.h", Sprintf("static uint64_t\t  base_mont[] = {%h, %h};\t\n// R = 2^{64 * words} mod p", Intseq(var_R, 2^64)[1], Intseq(var_R, 2^64)[2]));

var_R2 := var_R^2 mod var_p;
//print Intseq(var_R2, 2^64):Hex;
PrintFile("setup_FF.h", Sprintf("static uint64_t\tbase_mont_2[] = {%h, %h};\t\n// R^2 mod p", Intseq(var_R2, 2^64)[1], Intseq(var_R2, 2^64)[2]));

_, var_mu, _ := XGCD( (-var_p) mod var_r, var_r);
var_mu := var_mu mod var_r;
//print Intseq(var_mu, 2^64):Hex;
PrintFile("setup_FF.h", Sprintf("static uint64_t\t\t mu[] = {%h, %h};\t\n// You must update mu in Arith.S: mu * p = -1 mod 2^{64 * words}", Intseq(var_mu, 2^64)[1], Intseq(var_mu, 2^64)[2]));

PrintFile("setup_FF.h", Sprintf("static uint64_t\t      zeroM[] = {0x0000000000000000, 0x0000000000000000};"));
PrintFile("setup_FF.h", Sprintf("static uint64_t\t\tuno[] = {0x0000000000000001, 0x0000000000000000};\t\n// \"Normal\" One for exponentiation"));
PrintFile("setup_FF.h", Sprintf("static uint64_t\t       unoM[] = {%h, %h};\t\n// Montgomery's one", Intseq(var_R, 2^64)[1], Intseq(var_R, 2^64)[2]));
PrintFile("setup_FF.h", Sprintf("static uint64_t\t       negM[] = {%h, %h};\t\n// Montgomery's minus one", Intseq(neg_R, 2^64)[1], Intseq(neg_R, 2^64)[2]));

PrintFile("setup_FF.h", "");
PrintFile("setup_FF.h", Sprintf("static uint64_t constant_sqrt1[] = {%h, %h};\t\n// Used for square root computation: (p - 3) / 4", Intseq( (var_p - 3) div 4, 2^64)[1], Intseq( (var_p - 3) div 4, 2^64)[2]));
PrintFile("setup_FF.h", Sprintf("static uint64_t constant_sqrt2[] = {%h, %h};\t\n// Used for square root computation: (p - 1) / 2", Intseq( (var_p - 1) div 2, 2^64)[1], Intseq( (var_p - 1) div 2, 2^64)[2]));

PrintFile("setup_FF.h", "");
PrintFile("setup_FF.h", Sprintf("#define WORD_N %o", var_words));
PrintFile("setup_FF.h", Sprintf("#define Log2_E %o\n", Ceiling(Log(2, var_p) / 2.0)));
// 8-bits words required for the kernel size
//PrintFile("setup_FF.h", Sprintf("#define FixedBYTES %o", Ceiling( Log(2^8, var_p) / 4.0)));
PrintFile("setup_FF.h", Sprintf("#endif"));

// Arith.S
PrintFile("Arith.new", Sprintf("prime:") : Overwrite:=true);
PrintFile("Arith.new", Sprintf("\t .quad  %h", Intseq(var_p, 2^64)[1]));
PrintFile("Arith.new", Sprintf("\t .quad  %h", Intseq(var_p, 2^64)[2]));
PrintFile("Arith.new", Sprintf("mu:") );
PrintFile("Arith.new", Sprintf("\t .quad  %h",   Intseq(var_mu, 2^64)[1]));
PrintFile("Arith.new", Sprintf("\t .quad  %h\n", Intseq(var_mu, 2^64)[2]));
