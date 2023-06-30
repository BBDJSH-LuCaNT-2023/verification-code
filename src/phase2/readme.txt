The file "phase1.c" is the code for phase 1. 
You should compile the file in the same directory as the archt source code, 
and link it against mpfi, mpfr and gmp. (Always make sure that you have the latest version of mpfi.)
Here's an example to illustrate the usage:

./phase1 0 1.84923950291870468835054358789774149509195663778321840e2 2.5e-51 333

The first number is the parity (0 or 1). 
The second and third numbers are the eigenvalue and proven error bound for it, 
which you can take from the R line of each of the cc_* files. 
The fourth number is the desired final precision of the eigenvalue, in bits.

The following is an example to output result from screen to a file.

e.g.	./phase1 ..>output.txt
If your shell is bash then you similarly re-direct stderr by putting a 2 directly of the '>'
e.g.	./phase1 ..>output.txt 2>/dev/null




