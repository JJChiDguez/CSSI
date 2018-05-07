//===========================================================
//=====================Isogenies=============================
//===========================================================

//Weiestrass 2-Isogenies 
get_2_isog := function(P, Curve)
    
    Px := P[1];

    v := 3*Px^2 + Curve[1]; 
    w := v*Px;
    
    isog := [v, Px];
    Curve_res := [Curve[1] - 5*v, Curve[2] - 7*w];
    return Curve_res, isog ; //S + M + 3Mc + 3A
end function;


eval_2_isog := function(P, isog)
    x := P[1];
    y := P[2];
    v := isog[1];
    xQ := isog[2];
    
    c := (x - xQ)^-1;
    c2 := c^2;
    
    X := v * c  + x ;
    beta := - y * ( v * c2 - 1);
    
    return [X, beta, 1]; // 1I + 1S +  3M  + 2A + 1Ac
end function;

//Weiestrass 3-isogenies

get_3_isog := function(Q, Curve)
    
    a := Curve[1];
    b := Curve[2];
    
    Px := Q[1];
    Py := Q[2];
    
    gPx := 3 * Px^2  + a;
    gPy := - 2 * Py ;
    
    v := 2 * gPx ;
    u := gPy^2;
    
    w := u + Px * v; 
    gPxy := gPx * gPy; //make this here instead in evaluation 
    
    Isog := [v, u, gPxy, Px, Py];
    
    Curve_New := [a - 5*v,  b  - 7*w];
    
    return Curve_New , Isog;
    
end function;//2S + 2M + 5MC + 4A

eval_3_isog := function(P, Isog)

    v  :=Isog[1];
    u  :=Isog[2];
    gQxy := Isog[3];

    xQ  := Isog[4];
    yQ  := Isog[5];
    
    x := P[1];
    y := P[2];
    
    c := (x - xQ)^-1;
    c2 := c^2;
    
    X :=  v * c  + u * c2 ;
    
    beta := ( (2 * u  * y  *   c)  + ( v * ( y - yQ ) )  -  gQxy   ) * c2;

    return [X + x , y - beta, 1];

end function; // 1I + 1S + 6M +1MC + 4A

//===========================================================
//=====================Get Points============================
//===========================================================

All2Torsion := function(E)
    Rts := Roots(X^3 + E[1]*X + E[2]);
    return Rts[1], Rts[2], Rts[3];
end function;

NoDual2Torsion := function(E, gamma)
    gamma2 := gamma^2;
    Sq := Sqrt(gamma2 - 4 * (E[1] + gamma2) );
    Inv2 := (F_q!2)^-1;
    return Inv2 * (-gamma + Sq), Inv2 * (-gamma - Sq);
end function;


//===========================================================
//====================Isomorphism============================
//===========================================================

get_isomorphism := function(E0, E1)
    if E0[2] eq 0 then
        u_aux4 := E1[1]*(E0[1]^(-1));
        u_aux2 := Sqrt(F_q!u_aux4);
        u_aux := Sqrt(u_aux2);
    else 
        u_aux2 := (E0[1] * E1[2]) * ( (E0[2] * E1[1])^-1 );
        u_aux := Sqrt(u_aux2);
    end if;
    
    return [u_aux, u_aux2];
end function;//2Sq + I + M || 4M + I + Sq

eval_isomorphism := function(P, Isom);
    return [Isom[2]*P[1], Isom[1] * Isom[2] * P[2]];
end function; //3M


//===========================================================
//====================Dual Isogeny===========================
//===========================================================

//2Dual======================================================
Dual_preKer := function(E, gamma)
    gamma2 := gamma^2;
    Sq := Sqrt(gamma2 - 4 * (E[1] + gamma2) );
    Inv2 := (F_q!2)^-1;
    return Inv2 * (-gamma + Sq);
end function; // 1S + 

get_2_dual := function(phi, D, C)
    //Input:
    //  D and C and a
    //  2-isogeny phi := D -> C
    //Output:
    //  Isogeny phi_d := C -> D
    //  Such that phi_d(phi) = [2]: D -> D and phi(phi_D) = [2] : C -> C.
    
    //Computing the kernel of the Dual isogeny phi_d
    preK := Dual_preKer(D, phi[2]);
    K1 := eval_2_isog([preK, 0], phi);
    Ed, phi_pred := get_2_isog(K1, C);
    Isom := get_isomorphism(Ed, D);
    return [phi_pred, Isom];
end function;

eval_2_dual := function(P, phi_d)
    P_d := eval_2_isog(P, phi_d[1]);
    return eval_isomorphism(P_d, phi_d[2]);
end function; // 1I + 1S +  6M  + 3A 
//Total Cost : 


//===========================================================
//====================Scalar Mul.============================
//===========================================================

//Point Addition
Pt_ADD := function(P, Q)
    xP := P[1];
    yP := P[2];
    xQ := Q[1];
    yQ := Q[2];
    
    if [xP, yP] eq [0,1] then 
        return [xQ, yQ];
    end if;
    
    if [xQ, yQ] eq [0, 1] then
        return [xP, yP];
    end if;
    
    L_U := yQ - yP;
    L_D := (xQ - xP)^-1;
    
    L := L_U * L_D;
    xR := L^2 - xP - xQ;
    yR := L * (xP - xR) - yP;
    
    return [xR, yR];
end function;

//Given P computes [-]P
Pt_Neg := function(P)
    return [P[1], -P[2]];
end function;

//Classical Point Doubling
Pt_DBL := function(P, a)

    X := P[1];
    Y := P[2];
    
    if [X, Y] eq [0, 1] then
        return [0, 1];
    end if;
    if Y eq 0 then
        return [0, 1];
    end if;
    
    l1 := 3* (X^2) + a;
    l2 := (2* Y )^(-1);
    L := l1 * l2;
    X2 := (L^2) -  (2 * X);
    Y2 := L * (X - X2) - Y;
    return [X2, Y2];

end function; // 1I + 2S + 2M + 3MC + 4A

//Computes [2^m] P
Pt_MDBL := function(P, a, m)
    for i := 1 to m do
        P := Pt_DBL(P, a);
    end for;
    return P;
end function;

//Computes [3] P
Pt_TRPL := function(P, a)
    //Perform P + 2P simultaneously
    X := P[1];
    Y := P[2];
    
    A := 3 * (X^2) + a;     
    A2 := A^2;             
    
    C := 2*Y;                   
    C2 := C^2;              
    
    B := A2 - 3*X*C2;  

    L_U := -(A * B + C2^2 );     
    L_I := (C * B)^-1;             
    L_b := (L_U * L_I);             
    
    C_I := L_I * B;
    
    X3 := L_b^2  - (A2 * (C_I^2) - X);
    Y3 := L_b * (X - X3) - Y;
    
    return [X3, Y3];
end function; // 1I + 5S + 6M + 3MC + 6A

//Computes [3^m] * P
Pt_MTRPL := function(P, a, m)
    for i := 1 to m do
        P := Pt_TRPL(P, a);
    end for;
    return P;
end function;

//Double and add algorithm to compute [m]P
DBLADD := function(P, m, a)
    Q := [0, 1];
    N := P;
    M_bin := Intseq(m, 2);
    for bit in M_bin do
        if bit eq 1 then
            Q := Pt_ADD(Q, N);
        end if;
        N := Pt_DBL(N, a);
    end for;
    return Q;
end function;

//===========================================================
//====================== Other ==============================
//===========================================================

JInvariant_2 := function(E_Curve);
    Jconst := 1728;
    a := E_Curve[1];
    b := E_Curve[2];
    U := 4 * a^3;
    D := (U + 27 * b^2 ) ^-1;
    
    return Jconst * U *D;
end function;



//==========================================================
    
Chain_iso_2Eval := function(var_R, cvr_E, var_W0, var_W1, var_e)

    E_i := cvr_E;
    R_i := var_R;
    W_0i := var_W0;
    W_1i := var_W1;
    duals := [];
    
    for i := 0 to var_e - 2 do
        ker := Pt_MDBL(R_i, E_i[1], var_e - i - 1);  // (2^(var_e - i - 1))*R_i;
		
        E_i, phi_i   := get_2_isog(ker, E_i);
        R_i := eval_2_isog(R_i, phi_i);
		W_0i := eval_2_isog(W_0i, phi_i);
		W_1i := eval_2_isog(W_1i, phi_i);
        
    end for;
    
	 E_i, phi_i := get_2_isog(R_i, E_i);
	 W_0i := eval_2_isog(W_0i, phi_i);
	 W_1i := eval_2_isog(W_1i, phi_i);
	 
    return W_0i, W_1i;
end function;

//==========================================================
    
Chain_dual_2Eval := function(var_R, cvr_E, var_W0, var_W1, var_e)

    E_j := cvr_E;
    R_i := var_R;
    W_0i := var_W0;
    W_1i := var_W1;
    duals := [];
    
    for i := 0 to var_e - 2 do
        ker := Pt_MDBL(R_i, E_j[1], var_e - i - 1);  // (2^(var_e - i - 1))*R_i;
		
        E_i, phi_i   := get_2_isog(ker, E_j);
        duals[i + 1] := get_2_dual(phi_i, E_j, E_i); 
        
        R_i := eval_2_isog(R_i, phi_i);
        E_j := E_i;
    end for;
    
	 E_i, phi_i := get_2_isog(R_i, E_j);
        duals[var_e] := get_2_dual(phi_i, E_j, E_i);

	for i in [1 .. var_e] do
		W_0i := eval_2_dual(W_0i, duals[var_e + 1 - i]);
		W_1i := eval_2_dual(W_1i, duals[var_e + 1 - i]);
       end for;
    return W_0i, W_1i;
end function;

VelusFormula := function(C, E)
    //Initialization===========================================
    //a invariants of EllipticCurve([a1, a2, a3, a4, a6]) 
    ai := Eltseq(E);
    a1 := ai[1];
    a2 := ai[2];
    a3 := ai[3];
    a4 := ai[4];
    a6 := ai[5];
    
    C2 := [];
    R := [];
    S := [];
    v := 0;
    w := 0;
    vQ_Set :=[];
    uQ_Set :=[];
    xQ_Set :=[];
    yQ_Set := [];
    gQx_Set :=[];
    gQy_Set := [];
    
    //Step 2 =================================================
    // C <- C \ {O}
    O := E ! 0;
    if (O in C) then
        Exclude(~C, O);
    end if;
    
    if (#C gt 0) then
        //Step 3 and 4========================================
        for p in C do
            if (2*p eq O) then
                //Set cointaining all 2-torsion points in C
                Append(~C2, p);
            else
                //Set R <- C \ {C2}
                Append(~R, p);
            end if;
        end for;

        //Setp 5 ===============================================
        //Split R into two equals sets such that Rplus and Rminus
        //so that if a point P is in Rplus then -P is in Rminus
        Rplus := [];
        Rminus := [];
        
        for p in R do
            if (-p in Rminus) then
                Append(~Rplus, p);
            else
                Append(~Rminus, p);
            end if;
        end for;

        //Step 6 ===============================================
        S := Rplus cat C2;
        
        //Step 7 ==============================================
        for Q in S do
            xQ := Q[1];
            yQ := Q[2];
            Append(~xQ_Set, xQ);
            Append(~yQ_Set, yQ);
            gQx := 3*Q[1]^2 + 2*a2* xQ + a4 - a1*yQ;
            gQy := -2*yQ -a1*xQ - a3;
            
            Append(~gQx_Set, gQx);
            Append(~gQy_Set, gQy);
            
            if 2*Q eq O then
                vQ := gQx;
            else
                vQ := 2*gQx - a1*gQy;
            end if;
            Append(~vQ_Set, vQ);
            
            uQ := gQy^2;
            Append(~uQ_Set, uQ);
            v +:= vQ;
            w +:= uQ + xQ*vQ; 
        end for;
    end if;
    
    //Step 8 =================================================
    A1 := a1;
    A2 := a2;
    A3 := a3;
    A4 := a4 - 5*v;
    A6 := a6 - (a1^2 + 4*a2)*v - 7*w;
    
    Ec := EllipticCurve([A1, A2, A3, A4, A6]);
    
    //We represent an Isogeny as the sets of values that 
    //will be used for evaluation purposes
    
    Isog := [xQ_Set, yQ_Set, vQ_Set, uQ_Set, gQx_Set, gQy_Set];
    return Ec, Isog;
end function;


Eval_Velu := function(P, Isog, E)
    ai := Eltseq(E);
    a1 := ai[1];
    a2 := ai[2];
    a3 := ai[3];
    a4 := ai[4];
    a6 := ai[5];
    
    xQ  := Isog[1];
    yQ  := Isog[2];
    vQ  :=Isog[3];
    uQ  :=Isog[4];
    gQx := Isog[5];
    gQy := Isog[6];
    
    x := P[1];
    y := P[2];
    alpha := 0;
    beta := 0;
    //Step 9=============================================================
    for i := 1 to #xQ do
        c := (x - xQ[i])^-1;
        c2 := c^2;
        c3 := c2* c;
        alpha +:=  vQ[i] *c  + uQ[i] * c2 ;
        beta +:= uQ[i]*(2*y + a1*x + a3)*(c3) + vQ[i] * (a1*c + y - yQ[i])  * c2 + (a1*uQ[i] - gQx[i]*gQy[i]) * c2;
    end for;
    
    return [alpha + x, y - beta,1];
end function;
/* ------------------------------------------------------------------ */
extract_kernel := function(var_Ks)
	return var_Ks[#var_Ks], Prune(var_Ks);
end function;

/* ----------------------------------------------------------------- */

/* ----------------------------------------------------------------- *
 * The functions:                                                    *
 *                spliting_scalars_Left, and                         *
 *                spliting_scalars_Right.                            *
 * Computes the points:                                              *
 *                      R, [2^(e / 2)]R, [2^(e / 4)]([2^(e / 2)]R),  *
 *                      [2^(e / 8)]([2^(e / 4)]([2^(e / 2)]R)), ..., *
 *                      [2^(e - 1)]R.                                *
 * However, this two functions compute the two possible cases:       *
 * Floor(e / 2^i) or Ceiling(e / 2^i).                               *
 *                                                                   *
 * This is very important because if e is not a power of two, then   *
 * the number of scalar multiplications and isogeny evaluations in   *
 * the computation of a d^e - isogeny, can be unbalanced.            *
/ ------------------------------------------------------------------ */

spliting_scalars_Left := function(var_d, var_e, var_R)
	aux_R := var_R;
	aux_e := var_e; acc_e := var_e;
	seq_K :=[]; idx := 0;
	cost_SM := 0; 
	while(aux_e gt 1) do
		idx +:= 1;

		seq_K[idx] := <aux_R, acc_e>;
		
		aux_R := var_d^(aux_e div 2) * aux_R;
		cost_SM +:= (aux_e div 2);
		
		acc_e -:= (aux_e div 2);
		aux_e := aux_e - (aux_e div 2);
	end while;
	
	seq_K[idx + 1] := <aux_R, acc_e>;
	return seq_K, cost_SM;
end function;

/* ----------------------------------------------------------------- */ 
spliting_scalars_Right := function(var_d, var_e, var_R)
	aux_R := var_R;
	aux_e := var_e; acc_e := var_e;
	seq_K :=[]; idx := 0;
	cost_SM := 0; 
	while(aux_e gt 1) do
		idx +:= 1;
		
		seq_K[idx] := <aux_R, acc_e>;
		
		aux_R := var_d^(aux_e - (aux_e div 2)) * aux_R;
		cost_SM +:= (aux_e - (aux_e div 2));
		
		acc_e -:= (aux_e - (aux_e div 2));
		aux_e := (aux_e div 2);
	end while;
	seq_K[idx + 1] := <aux_R, acc_e>;
	return seq_K, cost_SM;
end function;

/* ------------------------------------------------------------------ */
d_pow_e_isogenous_curve := function(var_R, var_E, var_d, var_e)
    E_i := var_E;
    R_i := var_R;
    cost_IE_Left := 0;
    
    /* ----------------------------------------------------------------- *
     * Here, we are computing the Left half of the trajectory in Isogeny *
     * Triangle.                                                         *
     * So, at the end of this half we have computed cost_SM_Lelft scalar *
     * multsiplications, and cost_IE_Left Isogeny evaluations.           *
     * ----------------------------------------------------------------- */
    
    var_Kernels, cost_SM_Left := spliting_scalars_Left(var_d, var_e, var_R);
    while(var_Kernels[1,2] ne ((var_e div 2)) ) do
    	
    	var_kernel, var_Kernels := extract_kernel(var_Kernels);
    	// If we have a kernel with order greater that var_d, then we spĺit such kernel, 
    	// and we add such splits in the list of kernels.    	
    	if( var_kernel[2] ne 1 ) then
	    	aux_Kernels, aux_cost_SM_Left := spliting_scalars_Left(var_d, var_kernel[2], var_kernel[1]);
    		
    		var_Kernels := var_Kernels cat aux_Kernels;
    		cost_SM_Left +:= aux_cost_SM_Left;
    	// Otherwise, the kernel has order var_d, and we can compute the var_d-isogeny.	
    	else
		C := [];
		for j := 1 to var_d do
			C[j] := j * var_kernel[1];
		end for;
    		E_i, phi_i := VelusFormula(C, E_i);
    		
    		aux_Kernels := [];
	    	// In addition, we evaluate each element in the sequence of Kernels, and we decrease
	    	// its order by 1.
    		for var_i in [1 .. #var_Kernels] do
    			aux_Kernels[var_i] := <E_i ! Eval_Velu(var_Kernels[var_i, 1], phi_i, E_i), var_Kernels[var_i, 2] - 1>;
    			cost_IE_Left +:= 1;
    		end for;
    		
    		var_Kernels := aux_Kernels;
    	end if;
    end while;
    
    /* ------------------------------------------------------------------ *
     * Here, we are computing the Right half of the trajectory in Isogeny *
     * Triangle.                                                          *
     * So, at the end of this half we have computed cost_SM_Right scalar  *
     * multsiplications, and cost_IE_Right Isogeny evaluations.           *
     * Morever, we have the following:                                    *
     *                                                                    *
     * cost_IE_Left = cost_SM_Right, and cost_IE_Right = cost_SM_Left.    *
     * ------------------------------------------------------------------ */
    //print "HEY YOU!";
    var_Kernels, cost_SM_Right := spliting_scalars_Right(var_d, var_Kernels[1,2], var_Kernels[1,1]);
    cost_IE_Right := 0;
    while(IsEmpty(var_Kernels) eq false) do

    	var_kernel, var_Kernels := extract_kernel(var_Kernels);
    	// If we have a kernel with order greater that var_d, then we spĺit such kernel, 
    	// and we add such splits in the list of kernels.
    	if( var_kernel[2] ne 1 ) then
	    	aux_Kernels, aux_cost_SM_Right := spliting_scalars_Right(var_d, var_kernel[2], var_kernel[1]);
    		
    		var_Kernels := var_Kernels cat aux_Kernels;
    		cost_SM_Right +:= aux_cost_SM_Right;
    	// Otherwise, the kernel has order var_d, and we can compute the var_d-isogeny.
    	else
		C := [];
		for j := 1 to var_d do
			C[j] := j * var_kernel[1];
		end for;
    		E_i, phi_i := VelusFormula(C, E_i);
    		
    		aux_Kernels := [];
	    	// In addition, we evaluate each element in the sequence of Kernels, and we decrease
	    	// its order by 1.
    		for var_i in [1 .. #var_Kernels] do
    			aux_Kernels[var_i] := <E_i ! Eval_Velu(var_Kernels[var_i, 1], phi_i, E_i), var_Kernels[var_i, 2] - 1>;
    			cost_IE_Right +:= 1;
    		end for;
    		
    		var_Kernels := aux_Kernels;
    	end if;
    end while;
    
    //printf ">> #Scalar mults. :\t%o\n", cost_SM_Left + cost_SM_Right; 
    //printf ">> #Isogeny eval. :\t%o\n", cost_IE_Left + cost_IE_Right;
    return E_i;
end function;

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
double_eval_d_pow_e_isogeny := function(var_R, var_W0, var_W1, var_E, var_d, var_e)
    E_i := var_E;
    R_i := var_R;
    cost_IE_Left := 0;
    W_0i := var_W0;
    W_1i := var_W1;
    
    /* ----------------------------------------------------------------- *
     * Here, we are computing the Left half of the trajectory in Isogeny *
     * Triangle.                                                         *
     * So, at the end of this half we have computed cost_SM_Lelft scalar *
     * multsiplications, and cost_IE_Left Isogeny evaluations.           *
     * ----------------------------------------------------------------- */
    
    var_Kernels, cost_SM_Left := spliting_scalars_Left(var_d, var_e, var_R);
    while(var_Kernels[1,2] ne (var_e - (var_e div 2)) ) do
    	
    	var_kernel, var_Kernels := extract_kernel(var_Kernels);
    	// If we have a kernel with order greater that var_d, then we spĺit such kernel, 
    	// and we add such splits in the list of kernels.    	
    	if( var_kernel[2] ne 1 ) then
	    	aux_Kernels, aux_cost_SM_Left := spliting_scalars_Left(var_d, var_kernel[2], var_kernel[1]);
    		
    		var_Kernels := var_Kernels cat aux_Kernels;
    		cost_SM_Left +:= aux_cost_SM_Left;
    	// Otherwise, the kernel has order var_d, and we can compute the var_d-isogeny.	
    	else
		C := [];
		for j := 1 to var_d do
			C[j] := j * var_kernel[1];
		end for;
    		E_i, phi_i := VelusFormula(C, E_i);
    		
    		aux_Kernels := [];
	    	// In addition, we evaluate each element in the sequence of Kernels, and we decrease
	    	// its order by 1.
    		for var_i in [1 .. #var_Kernels] do
    			aux_Kernels[var_i] := <E_i ! Eval_Velu(var_Kernels[var_i, 1], phi_i, E_i), var_Kernels[var_i, 2] - 1>;
    			cost_IE_Left +:= 1;
    		end for;
    		
		W_0i := E_i ! Eval_Velu(W_0i, phi_i, E_i);
		W_1i := E_i ! Eval_Velu(W_1i, phi_i, E_i);    		
    		var_Kernels := aux_Kernels;
    	end if;
    end while;
    
    /* ------------------------------------------------------------------ *
     * Here, we are computing the Right half of the trajectory in Isogeny *
     * Triangle.                                                          *
     * So, at the end of this half we have computed cost_SM_Right scalar  *
     * multsiplications, and cost_IE_Right Isogeny evaluations.           *
     * Morever, we have the following:                                    *
     *                                                                    *
     * cost_IE_Left = cost_SM_Right, and cost_IE_Right = cost_SM_Left.    *
     * ------------------------------------------------------------------ */
    
    var_Kernels, cost_SM_Right := spliting_scalars_Right(var_d, var_Kernels[1,2], var_Kernels[1,1]);
    cost_IE_Right := 0;
    while(IsEmpty(var_Kernels) eq false) do

    	var_kernel, var_Kernels := extract_kernel(var_Kernels);
    	// If we have a kernel with order greater that var_d, then we spĺit such kernel, 
    	// and we add such splits in the list of kernels.
    	if( var_kernel[2] ne 1 ) then
	    	aux_Kernels, aux_cost_SM_Right := spliting_scalars_Right(var_d, var_kernel[2], var_kernel[1]);
    		
    		var_Kernels := var_Kernels cat aux_Kernels;
    		cost_SM_Right +:= aux_cost_SM_Right;
    	// Otherwise, the kernel has order var_d, and we can compute the var_d-isogeny.
    	else
		C := [];
		for j := 1 to var_d do
			C[j] := j * var_kernel[1];
		end for;
    		E_i, phi_i := VelusFormula(C, E_i);
    		
    		aux_Kernels := [];
	    	// In addition, we evaluate each element in the sequence of Kernels, and we decrease
	    	// its order by 1.
    		for var_i in [1 .. #var_Kernels] do
    			aux_Kernels[var_i] := <E_i ! Eval_Velu(var_Kernels[var_i, 1], phi_i, E_i), var_Kernels[var_i, 2] - 1>;
    			cost_IE_Right +:= 1;
    		end for;
		W_0i := E_i ! Eval_Velu(W_0i, phi_i, E_i);
		W_1i := E_i ! Eval_Velu(W_1i, phi_i, E_i);
    		var_Kernels := aux_Kernels;
    	end if;
    end while;
    
    //	printf ">> #Scalar mults. :\t%o\n", cost_SM_Left + cost_SM_Right; 
    //	printf ">> #Isogeny eval. :\t%o\n", cost_IE_Left + cost_IE_Right;
    return W_0i, W_1i;
end function;

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

/*
K_L, C_L := spliting_scalars_Left(2, e_a, R_a);

for i in [ 1 .. #K_L] do
	printf "(%o, %o),", FactoredOrder(K_L[i,1])[1,2], K_L[i, 2];
end for;
printf "\n";

K_R, C_R := spliting_scalars_Right(2, e_a, R_a);

for i in [ 1 .. #K_R] do
	printf "(%o, %o),", FactoredOrder(K_R[i,1])[1,2], K_R[i, 2];
end for;
printf "\n";
*/

/* ------------------------------------------------------------------ */
chain_d := function(var_R, var_E, var_d, var_e)
    E_i := var_E;
    R_i := var_R;
    cost_IE_Left := 0;
    
    /* ----------------------------------------------------------------- *
     * Here, we are computing the Left half of the trajectory in Isogeny *
     * Triangle.                                                         *
     * So, at the end of this half we have computed cost_SM_Lelft scalar *
     * multsiplications, and cost_IE_Left Isogeny evaluations.           *
     * ----------------------------------------------------------------- */
    
    var_Kernels, cost_SM_Left := spliting_scalars_Left(var_d, var_e, var_R);
    while(IsEmpty(var_Kernels) eq false) do
    	
    	var_kernel, var_Kernels := extract_kernel(var_Kernels);
    	// If we have a kernel with order greater that var_d, then we spĺit such kernel, 
    	// and we add such splits in the list of kernels.    	
    	if( var_kernel[2] ne 1 ) then
	    	aux_Kernels, aux_cost_SM_Left := spliting_scalars_Left(var_d, var_kernel[2], var_kernel[1]);
    		
    		var_Kernels := var_Kernels cat aux_Kernels;
    		cost_SM_Left +:= aux_cost_SM_Left;
    	// Otherwise, the kernel has order var_d, and we can compute the var_d-isogeny.	
    	else
		C := [];
		for j := 1 to var_d do
			C[j] := j * var_kernel[1];
		end for;
    		E_i, phi_i := VelusFormula(C, E_i);
    		
    		aux_Kernels := [];
	    	// In addition, we evaluate each element in the sequence of Kernels, and we decrease
	    	// its order by 1.
    		for var_i in [1 .. #var_Kernels] do
    			aux_Kernels[var_i] := <E_i ! Eval_Velu(var_Kernels[var_i, 1], phi_i, E_i), var_Kernels[var_i, 2] - 1>;
    			cost_IE_Left +:= 1;
    		end for;
    		var_Kernels := aux_Kernels;
    	end if;
    end while;

    
    //printf ">> #Scalar mults. :\t%o\n", cost_SM_Left ; 
    //printf ">> #Isogeny eval. :\t%o\n", cost_IE_Left ;
    return E_i;
end function;
