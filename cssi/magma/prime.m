clear;

e_2 := 32;
ds := PrimesUpTo(1000);
ds := [ ds[i+2] : i in [1 .. (#ds - 2)] ];

while(e_2 le 64) do
	e_3 := Floor(e_2 / Log(2,3));
	p_plus_1 := 2^e_2 * 3^e_3;
	
	for di in ds do
		if(IsPrime(p_plus_1 * di - 1)) then Log(2, p_plus_1 * di - 1), e_2, e_3, di; break; end if;
	end for;
	
	e_2 +:= 2;
end while;


ds := PrimesUpTo(1000);
ds := [ ds[i+2] : i in [1 .. (#ds - 2)] ];
e_2 := 70;
e_3 := Floor((128 - e_2) / Log(2,3)) - 5;
p_plus_1 := 2^e_2 * 3^e_3;
	for di in ds do
	if(IsPrime(p_plus_1 * di - 1)) then Log(2, p_plus_1 * di - 1), e_2, e_3, di; break; end if;
end for;

ds := PrimesUpTo(1000);
ds := [ ds[i+2] : i in [1 .. (#ds - 2)] ];
e_2 := 80;
e_3 := Floor((128 - e_2) / Log(2,3)) - 5;
p_plus_1 := 2^e_2 * 3^e_3;
	for di in ds do
	if(IsPrime(p_plus_1 * di - 1)) then Log(2, p_plus_1 * di - 1), e_2, e_3, di; break; end if;
end for;
