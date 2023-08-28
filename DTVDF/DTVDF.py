
# A implementation of decentralized trapdoor verifiable delay function DTVDF
#
# run as /path/to/sage -python DTVDF.py
#
# Warning: unaudited code.
#
# Copyright (c) 2023 Wei Wang <weiwangscsc@bupt.edu.com>

import time
import math
from math import factorial

from Crypto.Util.number import inverse, size
from ecdsa.numbertheory import is_prime


import random
import hashlib

time_multiplier = 1000


def myfact(n :int):
    assert (n>=0)
    if n<2:
        return 1
    else:
        return n*myfact(n-1)

param = {
    # RSA modulus length, in bits.
    'rsa_modulus_length_in_bits': 1024,
    # Number of DTVDF shares needed to obtain a DTVDF result y
    # This is the threshold ω in the paper.
    'number_parties_needed': 2,
    # Number of players are corrupted in the scheme.
    'number_parties_corrupted':0,
    # Number of players engaging in the scheme.
    'number_parties_total': 3,
    # The factorial of n in the paper.
    'delta': myfact(3),
}

# The current version of the script.
__version__ = '2023-08-24/001'


urandom = random.SystemRandom()

def gcd(p: int, q: int) -> int:
    #Returns the greatest common divisor of p and q
    while q != 0:
        (p, q) = (q, p % q)
    return p

def random(n):
    return urandom.randint(0, n-1)  # inclusive 0 and n-1

def pow_mod(base, exponent, modulo):
    return pow(int(base), int(exponent), int(modulo))

def mod_1(a, b):
    return a % b

def mul_mod(op1, op2, modulo):
    return mod_1(int(op1) * int(op2), modulo)

def is_sophie_germain(n):
    return is_prime(n) and is_prime((n-1)/2)




#用的这个函数
def prime_gen_change2(param):
    #((pub_key, priv_key), ph1) = rsa.newkeys(1024)

    #N=p*q,where p=2*p1+1 and q=2*q1+1

    #N is 1024 bits
    p1=22477664609411115811764585530871397671644918508656968363395279165227515181191964394952975584329516345811304346752842235325165208926862467646321831306556109391609653
    q1=1149893074682290650662789001443543062233708223909929288384093692406801223515128786013496681854535378889342595764296587783880602765548521211213903

    #N is 2048 bits
    #p1 =1520760909515346305367990380626704828881533480985774821905369161325277043102054668743225163476654279819535334065276616324731538452062419866484041372115070215434324443089900433295858890143497159495300049281919505770228155315952735922014513248420995120629979767712193468186148470013074244537431136331846041542665360920615796132281
    #q1=4771988445224992977022630428658388959833772111982240594367901235076277699109819766861998446649705389415039582251659830799678010407778398554237775917426726850336949056027527530936129355490368990646462784510381677531328910675713942044300535533364066919743846425705811615088054884502150835563

    p = int(p1 * 2 + 1)
    q = int(q1 * 2 + 1)

    #print("p的值为：", p)
    #print("q的值为：", q)

    Nsize = size(p * q)
    print("Nsize:",Nsize)

    return (p, q)

#用的这个函数
def key_gen_change(param, primes):
    (p, q) = primes

    p2 = (p - 1) // 2
    q2 = (q - 1) // 2

    #print("p2的值为：", p2)
    #print("q2的值为：", q2)

    # φ(N)=(p-1)*(q-1)=4m.
    # Here, φ(N) is an RSA modulus, and φ(N) is Euler function of N.
    m = p2 * q2

    sk_unshared = {
        'p': p,
        'q': q,
        'm': m,
    }
    pk = {
        'n': p * q,
    }
    return (sk_unshared, pk)


def evaluate_poly(poly, point, m):
    ret = 0
    for i in range(len(poly)):
        ret = ret + mod_1((poly[i] * mod_1(point ** i, m)), m)
    return mod_1(ret, m)


def split_shamir_change(secret, number_coeffs, number_shares, modulus):
    a = [0] * number_coeffs
    a[0] = secret

    for i in range(1, number_coeffs):
        a[i] = random(modulus)
    delta_ni = inverse(param['delta'], modulus)

    s = [0] * number_shares
    for i in range(number_shares):
        s[i] = int(evaluate_poly(a, i + 1, modulus))
        s[i] = pow_mod(int((s[i] * delta_ni)), 1, modulus)
    return s,a


def deal_change(param, sk_unshared, pk,e222):
    # Generate shares for the secret key by Shamir splitting
    # and shares of the verification key.
    n =pk['n']

    s,a = split_shamir_change(secret=e222,
                     number_coeffs=param['number_parties_needed'],
                     number_shares=param['number_parties_total'],
                     modulus=sk_unshared['m'])

    # verification keys
    v_pre = random(n)
    assert(gcd(v_pre, n) == 1)
    v = mul_mod(v_pre, v_pre, n)

    vs = [0] * param['number_parties_total']
    for i in range(len(vs)):
        vs[i] = pow_mod(v, s[i], n)

    cv = [0] * param['number_parties_needed']
    for i in range(len(cv)):
        cv[i] = pow_mod(v, a[i], n)

    sk_shared = {
        'v': v,
        's': s,
        'vs': vs,
    }
    return sk_shared




def DTIVDF_y_shares(param, pk, sk_shared, message):
    xi = [0] * param['number_parties_total']
    for i in range(param['number_parties_total']):
        exponent = 2 * sk_shared['s'][i]
        # print("param['delta']:", param['delta'])
        # print("sk_shared['s'][i]:", sk_shared['s'][i])
        # print("message:",message)
        # print("exponent:",exponent)
        xi[i] = pow_mod(message, exponent, pk['n'])
    return xi


def lagrange(S, i, j, delta):
    ret = delta
    for j_prime in S:
        if j_prime != j:
            ret = (ret * (i - j_prime)) / (j - j_prime)
    return ret



def reconstruct_DTVDF_shares(param, pk, sigshares, message):
    n = pk['n']

    delta = param['delta']

    # 计算论文中的w
    w = 1
    quorum = list(range(1, param['number_parties_needed'] + 1))
    for i in quorum:
        exponent = 2 * lagrange(quorum, 0, i, delta)
        part = pow_mod(sigshares[i - 1], exponent, n)
        w = mul_mod(w, part, n)

    return w


def hash_transcript(**transcript):
    hexdigest = hashlib.sha256(str(transcript).encode('utf-8')).hexdigest()
    return int(hexdigest, base=16)

def lift_message(message, delta, n):
    return pow_mod(message, 4*delta, n)



def construct_proofs(param, pk, sk_shared, message, sigshares):
    n = pk['n']
    v = sk_shared['v']
    L = param['number_parties_total']
    xt = lift_message(message, 1, n)
    proofs = [0] * L
    quorum = list(range(L))
    for i in quorum:
        r = random(n)
        c = hash_transcript(script_version=__version__,
                            param=param,
                            pk=pk,
                            party_index=i,
                            v=v,
                            xt=xt,
                            vi=sk_shared['vs'][i],
                            xi2=pow_mod(sigshares[i], 2, n),
                            vp=pow_mod(v, r, n),
                            xp=pow_mod(xt, r, n))


        z = int(sk_shared['s'][i])*c + r
        proofs[i] = (z, c)

    return proofs



def verify_proofs(param, pk, sk_shared, proofs, message, sigshares):
    n = pk['n']
    v = sk_shared['v']
    xt = lift_message(message,1, n)

    quorum = list(range(param['number_parties_total']))

    for i in quorum:
        their_z, their_c = proofs[i]

        vp1 = pow_mod(v, their_z, n)
        vp2 = pow_mod(sk_shared['vs'][i], -their_c, n)

        xp1 = pow_mod(xt, their_z, n)
        xp2 = pow_mod(sigshares[i], -2*their_c, n)

        our_c = hash_transcript(script_version=__version__,
                                param=param,
                                pk=pk,
                                party_index=i,
                                v=v,
                                xt=xt,
                                vi=sk_shared['vs'][i],
                                xi2=pow_mod(sigshares[i], 2, n),
                                vp=mul_mod(vp1, vp2, n),
                                xp=mul_mod(xp1, xp2, n))
        assert(our_c == their_c)




def validate_param(param):
    assert(param['number_parties_needed'] >=
           param['number_parties_corrupted']+1)
    assert((param['number_parties_total'] - param['number_parties_corrupted'])
           >= param['number_parties_needed'])
    param['delta'] = factorial(param['number_parties_total'])


def setup(param, pem_file=None):
    validate_param(param)

    (sk_unshared, pk) = key_gen_change(param, prime_gen_change2(param))


    return (sk_unshared, pk)

def random_message(pk):
    return random(pk['n'])



def DTVDF_Yshare_Gen(param, pk, sk_shared, message_to_sign):
    # Construct the DTVDF share
    y_shares = DTIVDF_y_shares(param, pk, sk_shared, message_to_sign)
    # Construct the proof of DTVDF shares
    proofs = construct_proofs(param, pk, sk_shared,
                              message_to_sign, y_shares)
    return (y_shares, proofs)



def test_DTVDF_algorithms():

    (sk_unshared, pk) = setup(param)

    DTVDF_x = 15619920774592561628351138998371642294622340518469892832433140464182509560910157

    # 1111111111111KshareGen------------------
    print("")
    print("(1) Start testing KshareGen-------------------------")
    for _ in range(1):
        start_t = time.time() * time_multiplier

        # Quickly calculate y, where ph1 is equal to 4m
        phi_N = 4 * sk_unshared['m']
        e = pow(2, T, phi_N)

        # Generate of shares of the private key
        sk_shared = deal_change(param, sk_unshared, pk, e)
        KShareGen_time = round(time.time() * time_multiplier - start_t)
        print("    - KShareGen_time:", (KShareGen_time), "ms")


    #222222222222YshareGen------------------
    print("")
    print("(2) Start testing YshareGen-------------------------")
    for _ in range(1):
        start_t = time.time() * time_multiplier

        # Generate shares of DTVDF result y and generate the corresponding proof
        (DTVDF_y_shares, proofs) = DTVDF_Yshare_Gen(param, pk, sk_shared, DTVDF_x)

        YShareGen_time = round(time.time() * time_multiplier - start_t)
        print("    - yShareGen_time:", (YShareGen_time), "ms")



    #33333333333333333YshareVerify------------------
    print("")
    print("(3) Start testing YshareVerify-------------------------")
    for _ in range(1):
        start_t = time.time() * time_multiplier

        # Verify the proof
        verify_proofs(param, pk, sk_shared, proofs, DTVDF_x, DTVDF_y_shares)
        YshareVerify_time = round(time.time() * time_multiplier - start_t)
        print("    - YshareVerify_time:", (YshareVerify_time), "ms")

    #44444444444444444444444combine------------------
    print("")
    print("(4) Start testing combine-------------------------")
    for _ in range(1):
        start_t = time.time() * time_multiplier

        # Combine the DTVDF result y
        DTVDF_result_y_recombined = reconstruct_DTVDF_shares(param,
                                                             pk,
                                                             DTVDF_y_shares,
                                                             DTVDF_x)

        print("The combined DTVDF results y:",DTVDF_result_y_recombined)
        combine_time = round((time.time() * time_multiplier - start_t), 2)
        print("    - combine_time:", (combine_time), "ms")



if __name__ == '__main__':
    # The time difficulty parameter T of DTVDF
    T = 0
    for i in range(1):
        T = 1000 * (i + 1)
        print("")
        print("setting the time difficulty parameter", T)
        test_DTVDF_algorithms()
        print("OK")

