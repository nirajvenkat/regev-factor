# Regev's Algorithm
*Submission for QRISE 2024 - Microsoft Challenge*

[Video Description](https://youtu.be/o7NXxKnswUo)

Performs quantum resource estimation for [Regev’s Algorithm](https://arxiv.org/abs/2308.06572) for prime factorization and discrete logarithms.

This repository contains some key primitives of Regev's Algorithm, implem:
-  Gaussian state preperation
    - [gaussian.qs](gaussian.qs): Simplest implementation that calls Q# library `PreparePureStateD()`
    - [gaussian_mot.qs](gaussian_mot.qs): Arbitrary state preparation by [Möttönen et al.](https://arxiv.org/abs/quant-ph/0407010)
    - [gaussian_gr.qs](gaussian_gr.qs): Uses [Grover-Rudolph](https://arxiv.org/abs/quant-ph/0208112) state preparation which is meant specifically for probability distributions 
- Quantum modular exponentiation
    - [binary_exp.qs](binary_exp.qs): Binary exponentiation used in the [original paper](https://arxiv.org/abs/2308.06572)
    - [fib_exp.qs](fib_exp.qs): Fibonacci exponentiation introduced as an optimization by [Ragavan et. al.](http://arxiv.org/abs/2310.00899v3)