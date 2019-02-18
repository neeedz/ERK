--- load PHCpack
---
loadPackage("PHCpack", Configuration=>{"path"=>"/Users/Nida/Documents/phcpack/", "PHCexe"=>"./phc","keep files"=>true})
---
--- To get a start system to compute the mixed volume of a chemical reaction network
--- input the system augmented by conservation laws and
--- specialize all (non-zero) parameter values k_i and total concentrations c_i as 1.


--- full ERK network
R=CC[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12];

f1 = x11*x7-2*x1; 
f2 = x7*x9+x1-2*x2; 
f3 = x10*x7-2*x3; 
f4 = x12*x8-2*x4;
f5 = x10*x8+x4-2*x5; 
f6 = x8*x9-2*x6;
f7 = x1+x2+x3+x7-1; 
f8 = x4+x5+x6+x8-1;
f9 = -x7*x9-x8*x9+x2+x6;
f10 = -x10*x7-x10*x8+x3+x5; 
f11 = -x11*x7+x1+x5+x6; 
f12 = x1+x2+x3+x4+x5+x6+x9+x10+x11+x12-1;

ERK = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12};

mixedVolume(ERK)

--- MV(ERK) = 7

---fully irreversible ERK network
R=CC[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12];

f1 = x1 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 - 1;
f2 = x2 + x7 + x8 + x11 - 1;
f3 = x3 + x4 + x5 + x6 - 1;
f4 = x3*x12 - x4;
f5 = x4 - 2*x5;
f6 = x3*x9 - x6;
f7 = x11 - 2*x7;
f8 = x2*x10 - x8;
f9 = -x3*x9 + x7;
f10 = -x2*x10 + x5;
f11 = x1*x2 - x11;
f12 = -x3*x12 + x7 + x8;

irrERK = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12};

mixedVolume(irrERK)

---partialirrerk where the only reverse reaction is kon<>0
R=CC[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12];

f1 = x1 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 - 1;
f2 = x2 + x7 + x8 + x11 - 1;
f3 = x3 + x4 + x5 + x6 - 1; 
f4 = x12*x3-x4; 
f5 = x4-2*x5; 
f6 = x3*x9-x6; 
f7 = x2*x9+x11-2*x7; 
f8 = x10*x2-x8; 
f9 = -x2*x9-x3*x9+x7; 
f10 = -x10*x2+x5;
f11 = x1*x2-x11;
f12 = -x12*x3+x7+x8;

partialirrERK = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12};

---MV(irrERK) = 3 ?

---reduced ERK network
R = CC[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10];

f1 = x1+x3+x4+x5+x6+x7+x9+x10-1; 
f2 = x2+x3+x4-1;
f3 = x1*x2-x3;
f4 = x3-2*x4;
f5 = x2*x7-x5*x8+x4;
f6 = -x6*x8+x4;
f7 = -x2*x7+x10;
f8 = x8+x9+x10-1; 
f9 = x5*x8-x9;
f10 = -2*x10+x9;

redERK = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10};

mixedVolume(redERK)

--- MV(redERK) = 3