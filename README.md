---------------------------------------------------------

CLASSICAL MEET IN THE MIDDLE (MITM) ATTACK APPLIED FOR SOLVING THE CSSI PROBLEM


This research project was realized by:

	- Adj Gora,
	- Cervantes-Vázquez Daniel,
	- Chi-Domínguez Jesús-Javier,
	- Menezes Alfred, and
	- Rodríguez-Henríquez Francisco.

---------------------------------------------------------

INTRODUCTION	

This is a C-code implementation which solves the CSSI 
problem, and it was created by Jesús Javier Chi Domínguez 
<chidoys@gmail.com>, <jjchi@computacion.cs.cinvestav.mx>.

Special thanks to José Eduardo Ochoa Jiménez for its 
modular arithmetic implementation.

The elliptic curve arithmetic was obtained from the 
formulaes that Daniel J. Bernstein and Tanja Lange have 
in their website:
	https://www.hyperelliptic.org/EFD/oldefd/jacobian.html

In addition, this C-code includes the implementation of the 
following Meet in the middle attack (MITM):

		-) Naive approach (i.e., the MITM-basic version)
		-) Deph-First-Search (MITM-DFS) approach, and
		-) Van Oorschot-Wiener approach (VW).
		
To be precise, the framework of the CSSI problem to be 
solved is the following.

PUBLIC PARAMETERS (the points S_A, T_A, S_B, and T_B are randomly chosen):

	- p is a prime of the form 2^e_a*3^e_b*f - 1 with 2|e_a,e_b,
	- E / F_{p^2} : y^2 = x^3 + Ax + B is a supersingular curve in Weierstrass form,
	- E_A is a random [2^e_a]-isogenous curve to E (in Weierstrass form), 
	- <P_A, Q_A> = E(F_{p^2})[2^e_a] and <S_A, T_A> = E_A(F_{p^2})[2^e_a],
	- E_B is a random [3^e_b]-isogenous curve to E (in Weierstrass form), 
	- <P_B, Q_B> = E(F_{p^2})[3^e_b] and <S_B, T_B> = E_B(F_{p^2})[2^e_a],
	- the curve arithmetic and isogeny functions were done in Jacobian coordinates.

PRIVATE PARAMETERS:

	- The isogeny \phi_A that connects E and E_A and the secret kernel of \phi_A,
	- The isogeny \phi_B that connects E and E_B and the secret kernel of \phi_B

SOLUTION:
	
	- The methods MITM-basic, MITM-DFS, and VW find two isomorphic curves C_0 and C_1 such that

	ALICE CASE:
		C_0 = E / <[2^{e_a/2}]([n]P_A + Q_A)> or C_0 = E / <[2^{e_a/2}](P_A + [m]Q_A)>, and 
		C_1 = E_B / <[2^{e_a/2}]([k]S_A + T_A)> or C_1 = E_A / <[2^{e_a/2}](S_A + [l]T_A)>.

	In addition, the kernel of \phi_A is computed by evaluating the points P_A and Q_A 
	under the composition of these two degree-[2^{e_a/2}] isogenies and the isomorphism.
	(E,P_A,Q_A) |---> (C_0, P_{A,0}, Q_{A,0}) |---> (C_1, P_{A,1}, Q_{A,1}) |---> (E'_A ~ E_A, P'_A, Q'_A)
	And then by finding m_A (or n_A) such that P'_A = [m_A]Q'_A (or Q'_A = [n_A]P'_A), which 
	satisifies <P_A - [m_A]Q_A> = ker (\phi_A) (or <-[n_A]P_A + Q_A> = ker (\phi_A)).
	
	BOB CASE (not implemented yet):
		C_0 = E / <[3^{e_b/2}]([n]P_B + Q_B)> or C_0 = E / <[3^{e_b/2}](P_B + [m]Q_B)>, and 
		C_1 = E_B / <[3^{e_b/2}]([k]S_B + T_B)> or C_1 = E_B / <[3^{e_b/2}](S_B + [l]T_B)>.
		
	In addition, the kernel of \phi_B is computed by evaluating the points P_B and Q_B 
	under the composition of these two degree-[3^{e_b/2}] isogenies and the isomorphism.
	(E,P_B,Q_B) |---> (C_0, P_{B,0}, Q_{B,0}) |---> (C_1, P_{B,1}, Q_{B,1}) |---> (E'_B ~ E_B, P'_B, Q'_B)
	And then by finding k_B (or l_B) such that P'_B = [k_B]Q'_B (or Q'_B = [l_B]P'_B), which 
	satisifies <P_B - [k_B]Q_B> = ker (\phi_B) (or <-[l_B]P_B + Q_B> = ker (\phi_B)).
			
---------------------------------------------------------

PREREQUISITES:

Any version of gcc, and openmp library must be installed.
In particular, this version requires the Magma Computational 
Algebra System for configuring the instances to be solved. 

In next versions the Magma requirement will be replaced by 
any version of python (wait for it!)

---------------------------------------------------------

INSTALLING:

To configure the headers and assambly files

	- setup_FF.h, 
	- setup_EC.h, and
	- Arith.S
by running the following statement:	sh cofing.sh

The command sh cofing.sh creates a random instances to be 
solved. However, if you want to change the prime number p 
you should modify the file ./magma/ea.eb.f, whose contain 
should be in the following format:

	var_A := 32;	var_B := 20;	var_f := 23;

Now, to build the executables you shoud run the following 
statement (you need to have installed magma): sh build.sh

Until this step, you're able for solving the random instances.

---------------------------------------------------------

RUNNING:

	./bin/naive NUMBER_OF_CORES ALICE_OR_BOB
	./bin/dfs NUMBER_OF_CORES ALICE_OR_BOB
	./bin/lambda NUMBER_OF_CORES ALICE_OR_BOB W BETA NUMBER_OF_FUNCTIONS_TO_USED

Here, ALICE_OR_BOB can take the values A or B (but the 
case B is not implemented yet). In addition, the input 
NUMBER_OF_FUNCTIONS_TO_USED is optional and 
determines how many different versions of the function 
f you want to reached (only for estimations).

STATISTICAL EXPERIMENTS (only Alice case):

	sh run_naive NUMBER_OF_CORES NUMBER_OF_EXPERIMENTS
	python naive_results.py
	Output:
		- Measured running time,
		- Measured Clock Cycles,
		
	sh run_dfs NUMBER_OF_CORES NUMBER_OF_EXPERIMENTS
	python dfs_results.py
	Output:
		- Measured running time,
		- Measured Clock Cycles,
	
	sh run_lambda NUMBER_OF_CORES W BETA NUMBER_OF_EXPERIMENTS NUMBER_OF_FUNCTIONS_TO_USED
	python lambda_results.py
	Output:
		- Measured running time,
		- Measured running time for all dist. points reached,
		- Measured running time for all collisions reached,
		- Number of functions,
		- Measured Clock Cycles,
		- Distinguished points reached,
		- Different distinguished points reached,
		- Different collisions reached, and
		- Total collisions reached.

