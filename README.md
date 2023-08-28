# DTVDF

We implement a decentralized trapdoor verifiable delay function DTVDF based on Wesolowski’s VDF scheme (EUROCRYPT'2019). It is a function f: X → Y that can resist parallel computing while allowing for easy and quick verification of calculation results. The rapid computation of DTVDF results relies on the collaboration of multiple participants. 
  
 # Four algorithms of DTVDF
  
We implement KShareGen, YShareGen, YShareVer, and Combine four algorithms. DTVDF is designed based on a trapdoor VDF, a variant of verifiable secret sharing over the integers, and a variant of the non-interactive zero-knowledge proof protocol.

 # References
 
[Efficient Verifiable Delay Functions.](https://eprint.iacr.org/2018/623.pdf) Wesolowski, 2018
