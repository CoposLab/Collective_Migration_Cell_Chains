C/C++ codes for simulating models with input files from "Configuration".

To compile the code, use:
  1. Dd doublet:
       g++-13 -Wl,-ld_classic doublet_Dd.cpp -o EXEFILENAME
     
  2. MDCK-like doublet:
       g++-13 -Wl,-ld_classic doublet_MDCK.cpp -o EXEFILENAME
     

To run simulations, use:
  1. Stepping Dd doublet:
      ./EXEFILENAME -a 1.0 -b 1.0 -c 1.0 -e 0.16 -f 0.008 -g 10.5 -i 1.04 -j 1.0 -k 1.0 -l 1.0 -m 0.16 -n 0.006 -p 10.5 -q 1.04&

  2. Gliding Dd doublet:
      ./EXEFILENAME -a 1.0 -b 1.0 -c 1.0 -e 0.16 -f 0.008 -g 8.0 -i 1.04 -j 1.0 -k 1.0 -l 1.0 -m 0.16 -n 0.006 -p 8.0 -q 1.04&

  3. MDCK-like doublet:
      ./EXEFILENAME -a 1.0 -b 1.0 -c 1.0 -e 0.16 -f 0.008 -g 10.5 -i 20.0 -j 1.0 -k 1.0 -l 1.0 -m 0.16 -n 0.006 -p 10.5 -q 20.0&

User-specified parameters:
  -a & -j: membrane resting tension, gamma, of the leader & trailer
  -b & -k: membrane elastic stiffness, k, of the leader & trailer
  -c & -l: adhesive bonds stiffness, ka, of the leader & trailer
  -e & -m: bond formation rate, kon, of the leader & trailer
  -g & -p: threshold rupture load, Fcrit, of the leader & trailer
  -i & -q: protrusive strength, rho_1, of the leader & trailer

Note:
  Please remember to recompile if intercellular junction parameters (R, r, A, a) are changed in the cpp code.


List of C/C++ codes:
1. To simulate Dd tandem pair migration, use:

    doublet_Dd.cpp
   
2. To simulate MDCK-like tandem pair migration, use:

    doublet_MDCK.cpp
   
3. To introduce additional data structure, use:

    structs_mem.h

