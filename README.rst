About
=====
The repository contains input and output files of spectroscopy calculations done with different multiplet programs. In each folder there is an IPython notebook that summarizes the results of the calculations.

Multiplet programs
==================
Quanty
------
Binaries are provided on the project's `download page <http://www.quanty.org/download>`_. The version tested was compiled on 25/11/2018.

TTMult
------
The source code of the programs that make up TTMult are distributed via `BitBucket <https://bitbucket.org/cjtitus/ttmult/overview>`_. The programs were compiled using the default Makefile. The last commit considered was `27dba3c <https://bitbucket.org/cjtitus/ttmult/commits/27dba3c105c0bd26f3a0e9947c02d75847fb4842>`_.

Calculations
============
+-----+------+----------+------------+------+-----------------------------+
|     | Ion  | Symmetry | Experiment | Edge | Additional Details          |
+=====+======+==========+============+======+=============================+
| 001 | Fe2+ | C3v      | XAS, XMCD  | L2,3 | Crystal field, B[111]       |
+-----+------+----------+------------+------+-----------------------------+
| 002 | Fe2+ | C3v      | XAS, XMCD  | L2,3 | Crystal field, B[001]       |
+-----+------+----------+------------+------+-----------------------------+
| 003 | Fe2+ | Oh       | XAS, XMCD  | L2,3 | Ligand field, B[001]        |
+-----+------+----------+------------+------+-----------------------------+
| 004 | Fe2+ | Td       | XAS        | K    | Crystal field               |
+-----+------+----------+------------+------+-----------------------------+
| 005 | Fe2+ | Td       | XAS        | K    | Crystal field, 3d-4p mixing |
+-----+------+----------+------------+------+-----------------------------+

Timings
=======
The timings where obtained on machine with an Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz processor, running Debian 8.

