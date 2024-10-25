from random import Random
import math

import numpy as np
#from sage.matrix.matrix_space import MatrixSpace
#from  sage.rings.finite_rings.integer_mod_ring import  IntegerModRing
from sage.all import *
from utils import is_it_ternary, get_norm, addCN
from numpy.polynomial import  Polynomial
from fpylll import CVP, IntegerMatrix, GSO
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll import  BKZ as BKZ_FPYLLL
"""
Implementation of NTRU key generation for initial variant and NTREncrypt variant


How to:
    Create a NTRUKeyGenerator object with suitable parameters,
    generate seed via newSeed() and then call getKey() or get Lattice

Examples:

for dihderal group: 
    python attack.py 14  -q=128 --verbose=True --group="dihedral" --h="[115, 42, 117, 108, 73, 3, 53, 29, 108, 34, 72, 5, 36, 101]"

"""


class NTRUKeyGenerator:

    """
    Create NTRUKeyGenerator object with parameters
        -n: the order of the group n = 4*(smalln^2)
        -q: the required modulo
        -seed: the seed to generate a random bitstring (In the case of bqtru, it is difficult to generate f,g as required so the seed could be tuple)
        -h: the public key, one can specify public key or keep it empty

    """
    def __init__(self, n, q=0, seed=None, h=None, option=0, empty_fset= True, weak_instance= True, guess=False, attack_type=0):


        self.n   = n
        self.q   = q
        self.p   = 3
        self.empty_fset = empty_fset
        self.weak_instance = weak_instance
        self.guess = guess ## if False means: run the program without guessing
        self.attack_type = attack_type ## 0 means key recovery attack and 1 means message recovery attack
        self.option = option ## The option of doing the attack which determines the way of building the lattices
        if h==None:
            self.h = None
        else:
            self.h = h
        self.cache = {}
        self.seed = seed  ### it should be a tuple for bqtru
        self.sample_iid = 8*n
        self.sample_fixed_type = 30*n
        self.sample_key_bits = 2*self.sample_fixed_type
        self.FF  = IntegerRing()
        self.FFp = IntegerModRing(self.p)
        self.FFq = IntegerModRing(self.q)
        if self.n % 4 != 0:
            raise ValueError("The order should be multiple of 4: 4*(smallprime)^2!!!")
        self.N = int(self.n / 4)  ### For BQTRU n = 4*N = 4*smalln^2
        self.smalln = int(sqrt(self.N))
        if (round(self.smalln)-self.smalln!=0):
            raise  ValueError("The order should be 4*(smallprime)^2 !!!")
        # if self.n%2!=0:
        #     raise ValueError("The order of dihedral group is even!!!")
        if seed == None:
            self.seed_f = randint(0, 2 ** 64)
            self.seed_g = (randint(0, 2**64), randint(0, 2**64), randint(0, 2**64), randint(0, 2**64)) ##The seed of g is four values
            self.seed = ( self.seed_f, self.seed_g)  ## for bqtru, the seed is (seed_for_f, seed_for_g)
        else:
            self.seed_f = seed[0]
            self.seed_g = seed[1]
        self.d = math.floor((self.N/7))
        self.w = self.FFq(self.get_w())
        # print("q: ", self.q)
        # print("w: ", self.w)
        self.winv = self.w.inverse_of_unit()
        self.ninv = self.FFq(self.smalln).inverse_of_unit()  ##n^{-1} mod q
        self.nth_roots_of_unity = self.get_nth_root_of_unity()
        self.inv_nth_roots_of_unity = self.get_inv_nth_root_of_unity()
        self.nth_roots_for_Rq = self.get_nth_root_for_Rq()
        self.inv_nth_roots_for_Rq = self.get_inv_nth_root_for_Rq()
        self.lagrange_basis = self.get_lagrange_basis_two_variable()
        self.s = None
        self.guessed_s = None ## The guessed value of s
        self.T = None
        self.sigma_monomial = None
        self.sigma_lagrange = None
        self.V_monomial = None   ###delete later
        self.V_lagrange = None  ### delete later
        self.hstar = None
        self.f  = None
        self.g  = None
        self.Fq = None
        self.Fp = None
        self.D = None
        self.counting = 0
        self.guess_arr = [0] * self.N ### the guessing
        self.guess_arr[0] = 1 ###initially the set contains only the point(0,0)
        self.count = 1 ##initially the count is one
        self.counter = 0

        self.message       = None
        self.ciphertext    = None
        self.message_prime = None
        self.ciphertext    = None
        self.r             = None

        self.seedrm        = None  ### seed to generate m and r
        self.seedr         = None  ### seed to generate r
        self.seedm         = None  ### seed to generate m
        # we are going to access f,g in the attack class just to pass them and print them in a file

    """
    Input: Integer s.
    Output: Random bit array of length s.
    """

    def get_seed(self):
        """
        return the seed. For bqtru: the returned seed is a pair
        """
        return self.seed

    """
    Input: - seed for which we are generating a randomstring
           - s: the random seed for which the string is being generated
    """
    def randomBitArray(self,s, seed):

        random = Random()
        random.seed(seed)
        return [random.randrange(2) for i in range(s)]

    """
    set the seed as a random bit array of length sample_key_bits.
    To be used as seed in getKey() and getLattice().
    """
    def newSeed(self):
        # self.seed = self.randomBitArray(self.sample_key_bits, self.seed)
        self.seed  = randint(0, 2 ** 64)
        return self.seed

    """
    The seed for bqtru is generated separately and differently for f,g 
    """
    def newSeed_for_bqtru(self):

        # seed_f = self.randomBitArray(self.sample_fixed_type, self.seed_f)
        # seed_g = self.randomBitArray(self.sample_fixed_type, self.seed_g)
        seed_f = randint(0, 2 ** 64)

        seed_g0 = randint(0, 2 ** 64)
        seed_g1 = randint(0, 2 ** 64)
        seed_g2 = randint(0, 2 ** 64)
        seed_g3 = randint(0, 2 ** 64)

        seed_g = (seed_g0, seed_g1, seed_g2, seed_g3)
        self.seed = (seed_f, seed_g)
        return self.seed

    def update_seed(self, f=True):
        """
        :param f: True if we want to update the seed f
        otherwise, we update the seed g
        :return:
        """
        if f:
            self.seed_f = randint(0, 2**64)
        else:
            g0_seed = randint(0, 2**64)
            g1_seed = randint(0, 2**64)
            g2_seed = randint(0, 2**64)
            g3_seed = randint(0, 2**64)
            self.seed_g = (g0_seed, g1_seed, g2_seed, g3_seed)


    def set_seed_r(self, seedr):
        """

        :param seedr: an int in the range(0, 2**64) to generate a random sequence r used in the encryption
        :set self.seedr = seedr
        """
        self.seedr = seedr

    def set_seed_m(self, seedm):
        """

        :param seedm: an int in the range(0, 2**64) to generate a random sequence m for the message
        :return:
        """
        self.seedm = seedm


    def get_seedr(self):
        """

        :return: seedr
        """

        return  self.seedr

    def get_seedm(self):
        """

        :return: seedm
        """

        return  self.seedm

    def get_first_element_after(self, position, element=1):
        """

        Input: - position: we are looking for the first one after this position
        Returns the index of the first one after this position
        """

        for i in range(position, self.N):
            if self.guess_arr[i]!=(1-element):
                return i


    def get_first_element_before(self, position, element=1):
        """"
        Input: -position: we are looking for the first one before this position
        Returns the index of the first one before this position
        """
        for i in range(position-1, -1, -1):
            if self.guess_arr[i]!=(1-element):
                return i
        return -1

    def get_guess_as_set(self):
        """
        The function converts the guess_arr into a set
        """
        if self.guess_arr == 'over':
            return 'over'
        s = set({})
        for i in range(self.N):
            if self.guess_arr[i] != 0:
                s.add(i)
        return s

    def guess_next(self):
        """
        Call the guess function according to the correct instance
        """
        if self.weak_instance:
            self.get_next_guess_weak()
        else:
            self.get_next_guess_no_weak()
        return
    def get_next_guess_no_weak(self):
        """
        guess the set from which the set T has been constructed.
        This function returns the guess for the non-weak instances where the guessed set can be any set of C(n^k k) k<=n
        """
        self.counting +=1

        if self.counter < self.N-1 and self.guess_arr[self.counter+1] == 0:
            # print("branch 1")
            self.guess_arr[self.counter] = 0
            self.counter = self.counter+1
            self.guess_arr[self.counter] = 1

        else:
            acc = 0
            # print("branch 2")
            pos = self.get_first_element_before(self.counter, 1)
            while pos==self.counter-1:
                self.guess_arr[self.counter] = 0
                self.counter = self.counter-1
                acc +=1
                pos = self.get_first_element_before(self.counter, 1)
            # print("pos in branch 2: ", pos)
            if pos ==-1 :
                if self.count<self.smalln:
                    self.count = self.count+1
                    self.guess_arr =[0]*self.N
                    for i in range(self.count):
                        self.guess_arr[i] = 1
                    self.counter = self.count-1

                else:
                    self.guess_arr ='over'

            else: ##pos !=-1
                self.guess_arr[self.counter] = 0
                # print("pos: ", pos)
                # print("counter: ", self.counter)
                # print("count: ", self.count)
                self.guess_arr[pos] = 0
                self.counter = pos+1

                self.guess_arr[self.counter] = 1
                for i in range(acc+1):

                    self.counter = self.counter+1
                    self.guess_arr[self.counter] = 1
                #self.counter = self.counter-1



    def get_next_guess_weak(self):
        """
        guess the set from which the set T has been constructed.
        This function returns the guess for the weak instances; the way BQTRU is designed
        """
        self.counting +=1

        if self.count ==1:
            self.counter =1
            self.guess_arr[self.counter] = 1 ##the second list
            self.count = 2
            return
        if self.counter < self.N-1 and self.guess_arr[self.counter+1] == 0:
            # print("branch 1")
            self.guess_arr[self.counter] = 0
            self.counter = self.counter+1
            self.guess_arr[self.counter] = 1

        else:
            acc = 0
            # print("branch 2")
            pos = self.get_first_element_before(self.counter, 1)
            while pos==self.counter-1:
                self.guess_arr[self.counter] = 0
                self.counter = self.counter-1
                acc +=1
                pos = self.get_first_element_before(self.counter, 1)
            # print("pos in branch 2: ", pos)
            if pos ==0 :
                if self.count<self.smalln:
                    self.count = self.count+1
                    self.guess_arr =[0]*self.N
                    for i in range(self.count):
                        self.guess_arr[i] = 1
                    self.counter = self.count-1

                else:
                    self.guess_arr ='over'

            else: ##pos !=0
                self.guess_arr[self.counter] = 0
                self.guess_arr[pos] = 0
                self.counter = pos+1

                self.guess_arr[self.counter] = 1
                for i in range(acc+1):

                    self.counter = self.counter+1
                    self.guess_arr[self.counter] = 1
                #self.counter = self.counter-1



    """
       Input:  - l a list
               - modulo: a modulo
       Output: - The list l% modulo
    """
    def Modulo(self, l, modulo):

        if modulo == None:
            return l
        else:
            return [i % modulo for i in l]

    def center_lift_form(self,f, q=None):
        """
        Centerlifting a vector f with respect to modulo q
        Input: f is a list
               q: a modulo
        Output: the centerlifting of the vector f with respect to q
        """
        if q == None:
            q = self.q

        t = f[:]
        for i in range(len(f)):
            t[i] = int(t[i])
            if t[i] > int(q / 2):
                t[i] = t[i] - q
        return t


    """
     Input: two elements representing two polynomials from the ring Z_mod[x]/(X^n-1)
     Output: their multiplication
    """
    def ZCn_multiply(self, element1, element2, mod=None):


        N = len(element1)
        multi_result = [0] * N
        for i in range(N):
            for j in range(N):

                multi_result[(i + j) % N] = (multi_result[(i + j) % N] + element1[i] * element2[
                    j])  # ai*aj*ri*rj = ai*aj(r_(i+j))

        self.Modulo(multi_result, mod)

        return multi_result

    def scalar_by_list(self, sc, li):
        """
        Input: -sc: a scalar value to be multiplied by the list
               -li: a list
        """
        output_list = [sc * ele for ele in li]

        return output_list

    def addCN(self, a, b):
        """
        Input: a,b two elements in C_N
        Output: adding the corresponding coefficients.
        """
        n = len(a)
        c = [0] * n
        for i in range(n):
            c[i] = a[i] + b[i]
        return c

    def substractCN(self, a, b):
        """
        Input: a,b two elements in C_N
        Output: adding the corresponding coefficients.
        """
        n = len(a)
        c = [0] * n
        for i in range(n):
            c[i] = a[i] - b[i]
        return c

    def bivariate_multiply(self, element1, element2):
        """
        Input:  - n: defines x^n-1
                - element1, element2: two elements in R= Z[x,y]/(x^n-1, y^n-1)

        element = f0(x) + yf1(x)+ .....y^{n-1}f_{n-1}(x)

        Output: - the multiplication result
        """
        n = self.smalln
        multi_result = [0]*(n**2)
        for i in range(n):
            for j in range(n):
                index = (i + j) % n
                index = index * n
                multi_result[index:index + n] = self.addCN(multi_result[index:index + n],
                                                      self.ZCn_multiply( element1[i * n:i * n + n], element2[
                                                                                                 j * n:j * n + n]))  # ai*aj*ri*rj = ai*aj(r_(i+j)
        return multi_result




    def multiplication_in_A(self, qin1, qin2, modulo=None):
        """
        Input: q1, q2: two elements in A written as: q1 = a0+a1*i+a2*j+a3*k
        q1 = a0+a1*i+a2*j+a3*k is represented as [[a0], [a1], [a2], [a3]] and ai in R or R_p or R_q
        depending on the modulo value R= Z[x,y]/(x^n-1, y^n-1)

        Output: multiplying the two elements according to the rules of the quaternion
        i^2 = j^2 = 1
        k^2 = -1
        ij = -ji = k
        The coefficients themselves are just bivariate polynomials
        """
        q1 = copy(qin1)
        q2 = copy(qin2)

        N = int(len(q1) / 4)  ##n^2
        n = int(sqrt(N))  ##n
        q1 = [q1[:N], q1[N:2 * N], q1[2 * N:3 * N], q1[3 * N:4 * N]]
        q2 = [q2[:N], q2[N:2 * N], q2[2 * N:3 * N], q2[3 * N:4 * N]]
        constant = self.bivariate_multiply( q1[0], q2[0])  ##constant part
        constant = self.addCN(constant, self.scalar_by_list(1, self.bivariate_multiply( q1[1], q2[1])))  ## i^2 =  1
        constant = self.addCN(constant, self.scalar_by_list(1, self.bivariate_multiply(q1[2], q2[2])))  ## j^2 =  1
        constant = self.addCN(constant, self.scalar_by_list(-1, self.bivariate_multiply(q1[3], q2[3])))  ## k^2 = -1

        ipart = self.bivariate_multiply( q1[0], q2[1])  ### 1*i
        ipart = self.addCN(ipart, self.bivariate_multiply(q1[1], q2[0]))  ### i*1
        ipart = self.addCN(ipart, self.scalar_by_list(-1, self.bivariate_multiply( q1[2], q2[3])))  ### j*k = -i
        ipart = self.addCN(ipart, self.scalar_by_list(1, self.bivariate_multiply(q1[3], q2[2])))  ### k*j =i

        jpart = self.bivariate_multiply( q1[0], q2[2])  #### 1*j
        jpart = self.addCN(jpart, self.scalar_by_list(1, self.bivariate_multiply( q1[1], q2[3])))  #### i*k =j
        jpart = self.addCN(jpart, self.bivariate_multiply( q1[2], q2[0]))  ### j*1 = j
        jpart = self.addCN(jpart, self.scalar_by_list(-1, self.bivariate_multiply( q1[3], q2[1])))  ### k*i = -j

        kpart = self.bivariate_multiply(q1[0], q2[3])  ##1*k
        kpart = self.addCN(kpart, self.bivariate_multiply( q1[1], q2[2]))  ##i*j=k
        kpart = self.addCN(kpart, self.scalar_by_list(-1, self.bivariate_multiply( q1[2], q2[1])))  ##j*i=-k
        kpart = self.addCN(kpart, self.bivariate_multiply( q1[3], q2[0]))  ##k*1 = k

        return self.Modulo(constant + ipart + jpart + kpart, modulo)



    """
           From https://ntru.org/f/ntru-20190330.pdf, Section 1.10.5.
           Input: A bit array b of length sample_fixed_type_bits.
           Output: A ternary polynomial with exactly d1 coefficients equal to 1 and d2  coefficients equal to −1.
    """
    def fixed_type(self, b, d1, d2):
        A = [0] * (self.N )
        v = [0] * (self.N )
        i = 0
        while i < d1:
            A[i] = 1
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1
        while i < d1+d2:
            A[i] = 2
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1

        while i < self.N:
            for j in range(30):
                A[i] += 2 ** (2 + j) * b[30 * i + j]
            i += 1
        A.sort()
        for i in range(self.N):
            v[i] = A[i] % 4
            if v[i] ==2:
                v[i] =-1

        return v



    def ternary(self, b):
        """
        Input: - b: A bit array b of length sample_iid_bits.
               - n: the order of the group

        Output: A ternary polynomial.

       """
        v = [0] * self.N

        for i in range(self.N):
            coeff_i = 0
            for j in range(8):
                coeff_i += 2 ^ j * b[8 * i + j]
            v[i] = coeff_i

        for i in range(self.N):
            v[i] = v[i] % 3
            if v[i] == 2:
                v[i] = -1
        return v
    """
           Input: arr: an array , n: an integer 
           Output: shifting to left the array by n positions
    """

    def shiftLbyn(self, arr, n=0):
        return arr[n::] + arr[:n:]

    """
    
    Matrix representation for an element of Z_qC_n (right circulant matrix)
    It's also an auxiliary matrix for matrix representation for an element in Z_qD_n
    Input: the first row that represents an element f,g, or h.
    FF: the space over it, the matrix to be constructed either IntegerModRing(3)
    or IntgerModRing(q)
    

    Output: the matrix that represents the group ring element for a cyclic group
    """

    def get_A(self, first_row, FF):

        n = len(first_row)
        a = first_row
        m = []
        for i in range(n):
            m.append(a)
            a = self.shiftLbyn(a, -1)

        MS2 = MatrixSpace(FF, n, n)
        B = MS2.matrix(m)

        return B


    def multiplication_rules(self, index1, index2):
        """
        Input: index1, index2 the indices of the elements to be multiplied
        Output: the index of the output of the multiplication
        """
        indices = [1, 2, 3]
        if (index1 == 0):
            index = index2
        elif (index2 == 0):
            index = index1
        elif (index1 == index2):
            index = 0
        else:
            indices.remove(index1)
            indices.remove(index2)
            index = indices[0]
        return index



    ### We are showing the Ecludein part of BQTRU as a twisted Group Ring

    def get_lambda_bqtru(self, index1, index2):
        """
        Input: -n: the order of the group
               -index1, index2: the indicies of the element for some group
        """
        lamda = 1
        #     print("index 1:", index1)
        #     print("index 2: ", index2)

        if (index1 == 2):  ###for j
            if (index2 == 1 or index2 == 3):  ##ji or j^k
                lamda = -1


        elif (index1 == 3):
            if (index2 == 1 or index2 == 3):  ##ki or k^2
                lamda = -1  #

        return lamda


    def get_lambda_bqtru_tilde(self, index1, index2):
        """
        Input: -index1, index2: the indices for the group
        returns 1 or -1 according to the positions of the group elements so that
        we can get the matrix of the element F^{tilde} = [F0 F1 F2 F3 // F1 F0 -F3 - F2// F2 F3 F0 F1// -F3 - F2 F1 F0]
        """
        lamda = 1
        if index1 == 1:
            if index2==2 or index2 ==3:
                return -1

        if index1 == 3:
            if index2==2 or index2==3:
                return -1

        return lamda


    """
       Matrix representation for an element of BQTRU
       Input: the first row that represents an element f,g, or h.
       - FF: the space over it, the matrix to be constructed either IntegerModRing(3)
       or IntgerModRing(q)
       - tilde: if 0 we return the group ring matrix otherwise we return the tilde matrix
       where tilde matrix is the matrix used in the message recovery attack
       Output: the matrix that represents the group ring element of BQTRU ring
       """

    def get_BQTRU_mat(self, first_row, FF=IntegerRing(), tilde=0):

        N = self.N ## for bqtru N = n^2
        alpha_mat = self.get_bivariate_mat(first_row[:N], FF)
        beta_mat = self.get_bivariate_mat(first_row[N:2 * N], FF)
        gamma_mat = self.get_bivariate_mat(first_row[2 * N:3 * N], FF)
        delta_mat = self.get_bivariate_mat(first_row[3 * N:4 * N], FF)

        element_list = [alpha_mat, beta_mat, gamma_mat, delta_mat]
        mat = []
        a = []
        for i in range(4):
            for j in range(4):
                index = self.multiplication_rules(i, j)
                if tilde ==0:
                    lamda = self.get_lambda_bqtru(i, index)
                else:
                    lamda = self.get_lambda_bqtru_tilde(i, index)
                element = element_list[index]
                t = np.zeros(([element.nrows(), element.ncols()]))
                for j in range(element.nrows()):
                    for k in range(element.ncols()):
                        t[j][k] = lamda * element[j][k]
                MS2 = MatrixSpace(FF, element.nrows(), element.ncols())
                a.append((MS2.matrix(t)))
        M = block_matrix(4, a)
        return M

    """
      Input: - element: the element from quaternion
             - tilde: default 0 means we return the bqtru corresponding matrix
             - tilde =1, we return the matrix of the tilde 
            
       Returns the matrix corresponding to the element
       
    """
    def element_to_matrix(self,element, FF, tilde=0):

        return self.get_BQTRU_mat(element, FF, tilde)


    #
    # """
    #     Returns tilde matrix that is used for message recovery attack
    # """
    # def get_BQTRU_mat_tilde(self, element, FF=IntegerRing()):
    #
    #     BQTRU_mat = self.get_BQTRU_mat(element, FF) ##get BQTRU mat
    #     for i in range(self.N, 2*self.N):
    #         for j in range(2*self.N, 4*self.N):
    #             BQTRU_mat[i][j] = -1*BQTRU_mat[i][j]   ###flipping F3, F2
    #
    #     for i in range(2*self.N, 3*self.N):
    #         for j in range(self.N, 2*self.N):
    #             BQTRU_mat[i][j] = -1*BQTRU_mat[i][j]  #flipping F3
    #             BQTRU_mat[i][j+2*self.N] = -1*BQTRU_mat[i][j]  #flipping F2
    #
    #     for i in range(3*self.N, 4*self.N):
    #         for j in range(self.N, 2*self.N):
    #             BQTRU_mat[i][j] = -1*BQTRU_mat[i][j]  #flipping F2
    #             BQTRU_mat[i][j+self.N] = -1*BQTRU_mat[i][j] #flipping F1
    #
    #
    #     return




    ##########    functions defined for bqtru ###############

    def get_w(self):
        """
        The function returns the principle root of unity
        """
        flag = True
        for i in range(2, self.q - 1):
            if ((i**self.smalln) % self.q == 1):
                for j in range(1, self.smalln):
                    if (i**j % self.q== 1):
                        flag = False
                        break
                if flag == True:
                    return i

    def get_nth_root_of_unity(self):
        """
        Input: -w: the principle root of unity
        Output: - the function return the nth roots of unity
        """

        nth_root_of_unity = []
        for i in range(self.smalln):
            nth_root_of_unity.append(self.w**i)
        return nth_root_of_unity


    def get_inv_nth_root_of_unity(self,):
        """
        Input:  -w: the principle root of unity
        Output: - the function returns the inverse of the nth roots of unity
        """

        inv_nth_root_of_unity = []
        for i in range(self.smalln):
            inv_nth_root_of_unity.append(self.winv**i)
        return inv_nth_root_of_unity

    def get_nth_root_for_Rq(self):
        """
        Input: -w: the principle root of unity
        Output: the pairs of the nth root of unity for Rq
        """
        w = self.w
        nth_roots_for_Rq = []
        for i in range(self.smalln):
            for j in range(self.smalln):
                nth_roots_for_Rq. append((w**i, w**j))

        return nth_roots_for_Rq

    def get_inv_nth_root_for_Rq(self):
        """
        Input: -winv: the inverse of the  principle root of unity
        Output: the pairs of the inv nth root of unity for Rq
        """
        winv = self.winv
        inv_nth_roots_for_Rq = []
        for i in range(self.smalln):
            for j in range(self.smalln):
                inv_nth_roots_for_Rq.append((winv ** i, winv ** j))

        return inv_nth_roots_for_Rq

    def evaluate_one_variable(self, coeff, point):
        """
        Input:  - coeff: a list represents the coefficients of polynomial
                - point: value in Z_q
        Output: - the evaluation at this point
        """
        s = 0
        for i in range(len(coeff)):
            # print("coeff[i]: ", coeff[i])
            # print("point: ", point)
            s += coeff[i] * (point**i)
        return s % self.q

    def evaluate_two_variables(self, coeff, point):
        """
        Input: - coeff: a list represents the coefficients of a polynomial of two variables
               - point: a point to evaluate the function at
        """
        y_dim = int(len(coeff) / self.smalln)
        y_coef = []
        for i in range(y_dim):
            # print(point[0])
            y_coef.append(self.evaluate_one_variable(coeff[self.smalln * i:self.smalln * i + self.smalln], point[0]))  ##append the evaluation
        return self.evaluate_one_variable(y_coef, point[1])


    def get_lagrange_coeff(self, coeff, nth_roots_for_Rq):

        """
        Input: - coeff: a list represents the coefficients of a polynomial of two variables
               - nth_roots_for_Rq: a pair of points {(a,b) in Z_q x Z_q} represents a pair of points from the roots of unity
        Output: the coefficients of lagrange basis
        At this step: the coefficients are just simple evaluation
        """
        lagrange_basis = []
        for point in nth_roots_for_Rq:
            lagrange_basis.append(self.evaluate_two_variables(coeff, point))
        return lagrange_basis


    #### Another way form standard coeff to lagrange coeff that mimincs NTT
    def column(self, matrix, i):
        """
        Auxiliary function:
        Input: - matrix: the matrix
               - i: the index of the column to extract
        Output: the column of the matrix to extract
        """
        return [row[i] for row in matrix]

    def apply_NTT(self, coeff, roots_of_unity):

        """
        Input: - coeff(single variable function): a list represents the coefficients of polynomial
               - roots of unity: values in Z_q
               - q: the modulo

        Output: NTT a list that mimics the NTT evaluation for a function
        """
        NTT = []
        for root in roots_of_unity:
            NTT.append(self.evaluate_one_variable(coeff, root))  ##to be replaced in the future by more efficient method for evaluation
        return NTT


    def get_eval_mat(self, coeff, roots_of_unity):

        """
          Input: - coeff(two variables function): a list represents the coefficients of polynomial
                 - point: value in Z_q
                 - q: the modulo

          Output: matrix of NTT evaluations where every row represents evaluation for a single variable function
          mat = [NTT(f0),
                 NTT(f1),
                 NTT(f2),
                 ..
                 ..
                 NTT(f{n-1})
                  ]
          """
        mat = []
        y_dim = int(len(coeff) / self.smalln)
        for i in range(y_dim):
            mat.append(self.apply_NTT(coeff[self.smalln * i:self.smalln * i + self.smalln], roots_of_unity))
        return mat


    def get_final_eval(self, eval_mat, roots_of_unity):
        """
         Input: - eval_mat: the evaluation matrix
        eval_mat = [NTT(f0),
           NTT(f1),
           NTT(f2),
           ..
           ..
           NTT(f{n-1})
            ]
        roots_of_unity: the roots of unity
        n: the number of coefficients in one function(defines x^n-1)
        q: the modulo
        """

        lagrange_coeff = []
        for i in range(len(eval_mat)):
            lagrange_coeff += self.apply_NTT(self.column(eval_mat, i), roots_of_unity)

        return lagrange_coeff

    ############## get standard coeff from lagranage basis  #######################
    def get_standard_coeff(self, lagrange_coeff, inv_nth_roots):
        """
        Input: - lagrange_coeff: a list represents the coefficients of polynomial of two variables with respect to lagrange basis
               - n: the length of the one variable
               - ninv: 1/n in the ring FFq
               - nth_roots: a list contains n points corresponding to the nth roots of unity
               - q: the modulo

        Output: -The coefficients with respect to standard basis {1, x, x^2,.. x^{n-1}.. y, yx, ..yx^{n-1}}
        """

        ##first get eval mat
        n = self.smalln
        ninv = self.ninv
        # print("ninv: ", ninv)
        mat = []
        for i in range(n):
            t = self.apply_NTT(lagrange_coeff[n * i:n * i + n], inv_nth_roots)
            # print("apply NTT", t)
            mat.append(self.scalar_by_list(ninv, t))
        # print("inside standard mat ", mat)
        mat_tran = np.array(mat)
        mat_tran = mat_tran.transpose()

        mat = []
        for i in range(n):
            mat.append(mat_tran[i])

            ## At this point, the matrix has [NTT(f0); NTT(f1); .....]
        standard_coef = []
        for i in range(n):
            t = self.apply_NTT(mat[i], inv_nth_roots)
            standard_coef += self.scalar_by_list(ninv, t)

        return standard_coef

    ########## Getting lagrange basis ##################

    def get_lagrange_basis_one_variable(self, nth_roots):
        """
        Input: - n : defines x^n-1
               - nth_roots: the nth_roots of unity
               - FFq: the ring
        Output: polynomials corrsponding to lagrange basis
        """
        n = self.smalln
        FFq = self.FFq
        basis = []
        for i in range(n):
            b = 1
            for j in range(n):
                if i != j:
                    b = b * Polynomial([-1 * nth_roots[j], FFq(1)]) * FFq((nth_roots[i] - nth_roots[j]) ** -1)
            basis.append(b.coef.tolist())

        return basis

    def get_lagrange_basis_two_variable(self):
        """
        Input: - n : defines x^n-1
               - nth_roots: the nth_roots of unity
               - FFq: the ring
        Output: polynomials corresponding to lagrange basis
        """
        nth_roots = self.nth_roots_of_unity
        n = self.smalln
        FFq = self.FFq
        basis1 = self.get_lagrange_basis_one_variable(nth_roots)
        # print("basis1", basis1)
        basis2 = self.get_lagrange_basis_one_variable(nth_roots)
        # print("basis2", basis2)
        basis = []
        for i in range(n):
            for j in range(n):
                b = []
                for k in range(n):
                    b = b + self.scalar_by_list(basis2[j][k], copy(basis1[i]))
                basis.append(b)
        return basis


    def get_norm_bqtru(self, f):
        """
         Input:  n that defines x^n-1
                 f = f0+f1*i+f2j+f3k in A
         Output: the norm {f0}^2-{f1}^2-{f2}^2-{f3}^2
         """

        # print("inside the get norm bqtru: ", f)

        n = self.smalln
        # print("f[0:n^2]", f[0:n**2])
        f0s = self.bivariate_multiply(f[0:n**2], f[0:n**2])
        f1s = self.bivariate_multiply(f[n**2:2 * n**2], f[n ** 2:2 * n ** 2])
        f2s = self.bivariate_multiply(f[2 * n ** 2:3 * n ** 2], f[2 * n ** 2:3 * n ** 2])
        f3s = self.bivariate_multiply(f[3 * n ** 2:4 * n ** 2], f[3 * n ** 2:4 * n ** 2])
        s = self.substractCN(f0s, f1s)
        s = self.substractCN(s, f2s)
        s = self.addCN(s, f3s)

        return s

    def Rq_element_by_A(self, NF, fconjeguate):
        """
        Input:- NF: an element in R_q
              - fconjeguate: an element in A written as: s0+ s1i +s2j +s3k
        """

        ns = len(NF)  ## n^2: an element in Rq
        n = int(sqrt(ns))  ## n = sqrt(n^2)
        f0 = fconjeguate[0:ns]
        f1 = fconjeguate[ns: 2 * ns]
        f2 = fconjeguate[2 * ns: 3 * ns]
        f3 = fconjeguate[3 * ns: 4 * ns]
        f0p = self.bivariate_multiply(NF, f0)
        f1p = self.bivariate_multiply(NF, f1)
        f2p = self.bivariate_multiply(NF, f2)
        f3p = self.bivariate_multiply(NF, f3)
        return f0p + f1p + f2p + f3p


    def sample_g_modified(self, g_seed, lens, max_trials=100, weak_instance=True):
        """
         Input:  - len: the minimum length of the set T on it G_i vanishes
                 - max_trails to try before calling the function again
                 - weak_instance: True means the scheme is generated as in the paper otherwise no constraints on g
         Output:
                 - an element g ternary and vanishes on the set

        """

        # print("inside sample g modified function!")
        # print("weak instance inside sample g modified", weak_instance)
        s3 = None ## for the case of the weak instances
        seed0 = g_seed[0]
        seed1 = g_seed[1]
        seed2 = g_seed[2]
        seed3 = g_seed[3]

        first_loop =0
        n = self.smalln
        bitstr = 30 * n ** 2  ##number of coefficients is n^2
        while True:
            if first_loop!=0:
                seed0 = randint(0, 2 ** 64)
            else:
                first_loop = first_loop+1
            # print(seed)
            b = self.randomBitArray( bitstr, seed0)

            g0 = self.fixed_type(b, self.d, self.d)


            s0 = self.get_T(g0)
            # print("s0 for g0:  ", s0)
            if weak_instance == False:
                print("breaking for g0")
                break  ##no need to check any constraint for the weak_instance
            if len(s0) > lens:
                break

        trails = 0
        first_loop = 0
        while True:
            if first_loop!=0:
                seed1 = randint(0, 2 ** 64)
            else:
                first_loop = first_loop+1
            b = self.randomBitArray(bitstr, seed1)
            g1 = self.fixed_type( b, self.d, self.d)
            # print("before: ", weak_instance)
            # print("here : weak", ( weak_instance is False))
            # print("here : weak 2", (weak_instance is True))


            s1 = self.get_T(g1)
            if len(s1) >= lens:
                trails = trails + 1
                # print("s1 has increased !")
            # print("s1 for g1:  ", s1)
            s1 = s1.intersection(s0)
            if weak_instance == False:
                # print("breaking for g1")
                break  # no further check is needed

            if len(s1) >= lens:
                break
            if trails > max_trials:
                # print("calling sampling g function again!!")
                self.update_seed(False)
                g_seed = self.seed_g
                return self.sample_g_modified(g_seed, lens, max_trials)

        trails = 0
        first_loop = 0
        while True:
            if first_loop !=0:
                seed2 = randint(0, 2 ** 64)
            else:
                first_loop = first_loop+1

            b = self.randomBitArray( bitstr, seed2)
            g2 = self.fixed_type(b, self.d, self.d)


            s2 = self.get_T(g2)
            # print("s2 for g2: ", s2)
            if len(s2) >= lens:
                trails = trails + 1
                # print("trails have increased !")
            s2 = s2.intersection(s1)
            if weak_instance == False:
                # print("breaking for g2")
                break
            if len(s2) >= lens:
                break
            if trails > max_trials:
                # print("calling sample g function again!!")
                self.update_seed(False)
                g_seed = self.seed_g
                return self.sample_g_modified(g_seed, lens, max_trials)

        trials = 0
        while True:
            if first_loop!=0:
                seed3 = randint(0, 2 ** 64)
            else:
                first_loop= first_loop+1

            b = self.randomBitArray(bitstr, seed3)
            g3 = self.fixed_type( b, self.d, self.d)

            # print("s3 for g3: ", s3)
            s3 = self.get_T(g3)
            s3 = s3.intersection(s2)

            if weak_instance==False:
                # print("breaking for g3")
                break ##no further check is needed

            if len(s3) >= lens:
                trails = trails + 1

            if len(s3) >= lens:
                break
            if trails > max_trials:
                # print("calling sample g function again!!")
                self.update_seed(False)
                g_seed = self.seed_g
                return self.sample_g_modified(g_seed, lens, max_trials)

        self.seed_g = (seed0, seed1, seed2, seed3)  ###during the running of the program g_seed is being updated so we
                                                    ### we save the updated g_seed so we can return it
        return s3, g0 + g1 + g2 + g3

    def get_T(self, G):
        """
        Input:  - G: ternary element sampled as: g0+g1i+g2j+g3k where gi in T(d,d)
                - nth_roots_for_Rq: i.e. E = {all the points (a,b) in E}
                - q: the modulo

        Output: - T: the indices of the points for them G vanishes
        """

        # print("G inside get T ", G )
        nth_roots_for_Rq = self.nth_roots_for_Rq
        n = self.smalln
        s1 = set({})

        for i in range(n**2):
            # print("The following point is being evaluated: ", nth_roots_for_Rq[i])
            if self.evaluate_two_variables(G[0:n**2], nth_roots_for_Rq[i]) == 0:
                s1.add(i)
        # print("inside get_T: ", s1)
        #     print(s)
        if(len(G) == n**2):
            return s1      ## For the case of NF.

        s2 = set({})
        for i in range(n**2):
            if self.evaluate_two_variables(G[n**2:2 * n**2], nth_roots_for_Rq[i]) == 0:
                s2.add(i)
        # print("inside get_T: ", s2)

        s3 = set({})
        for i in range(n**2):
            if self.evaluate_two_variables(G[2 * n**2:3 * n**2], nth_roots_for_Rq[i]) == 0:
                 s3.add(i)
        # print("inside get_T: ", s3)

        s4 = set({})
        for i in range(n**2):
            if self.evaluate_two_variables(G[3 * n**2:4 * n**2],  nth_roots_for_Rq[i]) == 0:
                s4.add(i)
        # print("inside get_T: ", s4)

        s = s1.intersection(s2)
        s = s.intersection(s3)
        s = s.intersection(s4)
        return s


    def get_T_vectors(self):
        """
        Input: - lagrange_basis: the basis that defines the lagrange basis
               - the indices of T
        Output: - The T vectors of lagrange basis at the indices of s
        """
        lagrange_basis = self.lagrange_basis
        T = []
        for i in self.s:
            T.append(lagrange_basis[i])
        return T

    def get_geussed_T_vectors(self):
        """
        Input: - lagrange_basis: the basis that defines the lagrange basis
               - the indices of T
        Output: - The T vectors of lagrange basis at the indices of s
        """

        lagrange_basis = self.lagrange_basis
        T = []
        for i in self.guessed_s:
            T.append(lagrange_basis[i])
        return T


    # orginal_list = [1,x,x^2,...x^{n-1},  xy,           x^2y, ...x^{n-1}y, xy^2, ...x^n-1y]
    # inverse_list = [1,x^n-1, ..x,     , x^{n-1}y^{n-1},x^{n-2}y^{n-1} ]
    def get_inverse_bivariate(self, index):
        """
        Input: - index: the index of the group element (indexing start from zero)
               - n: is 4
        Output: the index of the inverse according to the group
        """

        n = self.N
        nprime = self.smalln
        if index == 0:
            return index
        else:
            ypos = int(index / nprime)
            xpos = index % nprime
            yinv = (nprime - ypos) % nprime
            xinv = (nprime - xpos) % nprime

            return (yinv * nprime + xinv) % n

    def multiplication_table_bivariate(self, index1, index2):
        """
        Input:  -n=2N: the order of the dihedral group
                -index1, index2: the indicies of the element for cyclic group
        Output: the index of the output when multiplying the group elements at positions index1, index2
        """
        n = self.N
        nprime = self.smalln
        ypos1 = int(index1 / nprime)
        xpos1 = index1%nprime

        ypos2 = int(index2 / nprime)
        xpos2 = index2 % nprime

        ypos = (ypos1 + ypos2) % nprime
        xpos = (xpos1 + xpos2) % nprime

        return ypos * nprime + xpos

    def get_bivariate_mat(self, element, FF=IntegerRing()):
        """
        Input:  -element: an element in R = Z[x,y]/(x^n-1, y^n-1)
        Output: -the matrix representation of the element in R

        """
        n = len(element)
        mat = []
        a = [0] * n
        for i in range(n):
            inv = self.get_inverse_bivariate(i)
            for j in range(n):
                index = self.multiplication_table_bivariate(inv, j)
                a[j] = element[index]
            # print(a)
            mat.append(copy(a))
        ###return matrix with entries from IntegerRing
        MS2 = MatrixSpace(FF, n, n)
        A = MS2.matrix(mat)
        return A

    def get_conjugate(self, s):
        """
            Input s: - an element in A
            Output sbar: the conjegute of s calculated as: s_0 - s_1i - s_2 j -s_3 k
        """
        length = int(len(s) / 4)
        return s[0:length] + self.scalar_by_list(-1, s[length:2 * length]) + self.scalar_by_list(-1, s[length * 2:length * 3]) + self.scalar_by_list(
            -1, s[length * 3:length * 4])




    def sample_fg_bqtru(self, f_seed, g_seed, f_set_empty=True, weak_instance=True):
        """
        Input: - f_seed, g_seed: two numbers in the range
               - f_set_empty: if True means fset doesn't need to vanish at any point
               - weak_instance: if True means the set T is defined where g vanishes as in BQTRU
        """
        # print("weak instance inside sample fg function: ", weak_instance)
        if weak_instance==False:

            f_set_empty=False ###T will be the f_set
            lens = 0 ##no constraint on the length of gset
        else:
            if self.smalln==5:
                lens = 1     ## for 5, so that decryption can work
            else:
                lens = 2     ## for other instances, len>=2 can work
        max_trails = 100
        s, g = self.sample_g_modified(g_seed, lens, max_trails, weak_instance)
        print("s after sampling g: ", s)
        # print("g seed after successfully sample g: ", g_seed)
        # g_bits = self.randomBitArray(self.sample_fixed_type, g_seed)
        # seqlen = int(self.sample_fixed_type / 4)
        # # g_bits = self.randomBitArray(self.sample_iid, g_seed)
        # # seqlen = int(self.sample_iid / 4)
        #
        #
        # ### sample G
        # g0 = self.fixed_type(g_bits[0:seqlen], self.d, self.d)
        # g1 = self.fixed_type(g_bits[seqlen:2 * seqlen], self.d, self.d)
        # g2 = self.fixed_type(g_bits[2 * seqlen:3 * seqlen], self.d, self.d)
        # g3 = self.fixed_type(g_bits[3 * seqlen:4 * seqlen], self.d, self.d)
        #
        # # g0 = self.ternary(g_bits[0:seqlen])
        # # g1 = self.ternary(g_bits[seqlen:2 * seqlen])
        # # g2 = self.ternary(g_bits[2 * seqlen:3 * seqlen])
        # # g3 = self.ternary(g_bits[3 * seqlen:4 * seqlen])
        #
        # g = g0 + g1 + g2 + g3
        # s = self.get_T(g)
        # print("s: ", s)
        # while (len(s)<=2 or len(s)>self.smalln):
        #     self.update_seed(False)  ##update g_seed
        #     print("updating seed g ")
        #     g_seed = self.seed_g
        #     g_bits = self.randomBitArray(self.sample_fixed_type, g_seed)
        #     # g_bits = self.randomBitArray(self.sample_iid, g_seed)
        #     ### sample G
        #     g0 = self.fixed_type(g_bits[0:seqlen], self.d, self.d)
        #     g1 = self.fixed_type(g_bits[seqlen:2 * seqlen], self.d, self.d)
        #     g2 = self.fixed_type(g_bits[2 * seqlen:3 * seqlen], self.d, self.d)
        #     g3 = self.fixed_type(g_bits[3 * seqlen:4 * seqlen], self.d, self.d)
        #     # g0 = self.ternary(g_bits[0:seqlen])
        #     # g1 = self.ternary(g_bits[seqlen:2 * seqlen])
        #     # g2 = self.ternary(g_bits[2 * seqlen:3 * seqlen])
        #     # g3 = self.ternary(g_bits[3 * seqlen:4 * seqlen])
        #
        #     g = g0 + g1 + g2 + g3
        #     s = self.get_T(g)

        # print("s for G: ", s)
        seqlen = int(self.sample_fixed_type / 4)
        while(True):
            f_bits = self.randomBitArray(self.sample_fixed_type, f_seed)

            ### sample F
            f0 = self.fixed_type(f_bits[0:seqlen], self.d+1, self.d)
            f1 = self.fixed_type(f_bits[seqlen:2 * seqlen], self.d, self.d)
            f2 = self.fixed_type(f_bits[2 * seqlen:3 * seqlen], self.d, self.d)
            f3 = self.fixed_type(f_bits[3 * seqlen:4 * seqlen], self.d, self.d)
            f = f0 + f1 + f2 + f3
            NF = self.get_norm_bqtru(f)
            fset = self.get_T(NF)

            if weak_instance == False:
                fsetsize = 3 #minimum size of the fset
            else:
                if f_set_empty==False:
                    fsetsize =1 ## weak instance but fset is not empty

            if f_set_empty==True:

                condition = len(fset)==0   ##sampling f with empty set
            elif weak_instance==False:
                condition = len(fset)<self.smalln ##in the case of the weak instance fset should be smaller than small n
            else:  ##fset=False and weak_instance=True
                condition = (fset.issubset(s)) and (len(fset) >= fsetsize)
            if (condition): ### and len(fset)>=1
                print("The condition of the set cardinality is satisfied !")
                mat = self.get_bivariate_mat(NF, IntegerModRing(self.p))
                if mat.is_invertible():
                    Fp = mat.inverse()[0]
                    print("The invertibility condition is satisfied !")
                    Fp = self.Rq_element_by_A(Fp, self.get_conjugate(f))
                    break
            self.update_seed(True)  ##update seed_f
            print("updating seed f")
            f_seed = self.seed_f
        if self.weak_instance==False:
            print("fset: ", fset)
            lower_bound = len(fset)
            if lower_bound==0:
                lower_bound =1
            s = set({})
            s = s.union(fset)
            Tsize = randint(lower_bound, self.smalln)
            print("Tsize: ", Tsize)
            ## generate random Tsize to be the size of T
            i =0
            while i<Tsize:
                index = randint(0, self.smalln**2-1)
                print("s: ", s)
                print("index: ", index)
                if s.issuperset({index}):
                    continue
                else:
                    s.add(index)
                    print("s after add: ", s)
                    i = i+1
        # print("gset: ", s)
        # print("G: ", g)
        print("s: ", s)
        print("fset: ", fset)
        # print("f: ", f)
        self.s = s
        Fq  = self.get_inverse_in_Aq(f, self.nth_roots_for_Rq, self.inv_nth_roots_of_unity)
        # print("check in the function: ", self.multiplication_in_A(f,Fq, self.q))
        ##### private values are not supposed to be saved, we save them for checking purposes ########
        self.f = f
        self.g = g
        self.Fq = Fq
        self.Fp = Fp
        ###############################################################################################
        #self.s = self.get_T(self.g)  ## The set that mentions the positions that defines sigma
        self.T = self.get_T_vectors()  ## The set of the points that defines sigma
        self.sigma_monomial, self.sigma_lagrange = self.get_sigma()
        self.seed_f = f_seed
        self.seed = (self.seed_f, self.seed_g)
        print("Final updated seed in sample fg bqtru: ", self.seed)
        return (f,g, Fp, Fq)  ## To be compatible with the other f, g function




    def get_sigma(self):

        """
        Input:  - n: defines the cryptosystem
                - qprime: the modulo
                - s: the indices of T
                - lagrange basis
        Output: -sigma written with respect to monomial basis as: sum{q_i*lambda(a,b)(x,y)} for (a,b) in T
                -sigma with respect to lagrange basis
        """
        s = self.s
        print("inside get_sigma:")
        print("self.s: ", self.s)
        lagrange_basis = self.lagrange_basis
        ns = self.smalln**2  ##The length of the vector
        t = [0] * ns
        for i in s:
            t[i] = randint(1, self.q - 1)  ##generate random nonzero qi

        monomial_basis = np.matmul(t, lagrange_basis)

        return monomial_basis, t

    def get_coeff_wise_inv(self, l):
        """
        Input: - l a list to be inverted
               - qprime: a prime value of q

        """
        s = []
        for i in range(len(l)):
            if l[i] == 0:
                s.append(0)
            else:
                s.append(self.FFq(l[i]).inverse_of_unit())

        return s

    def get_T_in_monomial_basis(self, s, lagrange_basis):
        """
        Input:  - n: defines the cryptosystm
                - s: the indices of T
                - lagrange basis
        Output: -T written with respect to monomial basis
        """

        ns = self.smalln ** 2  ##The length of the vector
        t = [0] * ns
        for i in s:
            t[i] = 1  ##put ones in the positions where s is not zero
        monomial_basis = np.matmul(t, lagrange_basis)
        return monomial_basis

    def mod_Q(self, f0, T):
        """
        Input: f0: an element in R_q represented with respect to lagrange coefficients
               T: the indices of T in E
        Output: f0 mod Q: which makes f0 largrange coefficients to be zeros at the indices of T
        """
        for i in T:
            f0[i] = 0
        return f0


    def mod_J(self, f, T):
        """
        Input: - f an element in A can be written as f0 + f1*i + f2*j +f3*k represented with respect to lagrange coefficients
               - T: the indices of T in E
               - n: the parameter that defines R
        Output: f mod J
        """
        n = self.smalln
        f0 = self.mod_Q(f[0:n ** 2], T)
        f1 = self.mod_Q(f[n ** 2: 2 * n ** 2], T)
        f2 = self.mod_Q(f[2 * n ** 2: 3 * n ** 2], T)
        f3 = self.mod_Q(f[3 * n ** 2: 4 * n ** 2], T)

        return f0 + f1 + f2 + f3

    def get_lagrange_coeff_for_A(self, A, nth_roots_for_Rq):
        """
        Input:  - A an element in the quaternion can be written as A0+A1*i+A2*j+A3*k
                - n, qprime, nth_roots_for_Rq: the parameters that define the cryptosystem
        Output: - The lagrange coefficients of  A = (A0, A1, A2, A3)
        """
        n = self.smalln
        A0 = self.get_lagrange_coeff(A[0:n ** 2], nth_roots_for_Rq)
        A1 = self.get_lagrange_coeff(A[n ** 2:2 * n ** 2],  nth_roots_for_Rq)
        A2 = self.get_lagrange_coeff(A[2 * n ** 2:3 * n ** 2],  nth_roots_for_Rq)
        A3 = self.get_lagrange_coeff(A[3 * n ** 2:4 * n ** 2], nth_roots_for_Rq)

        return A0 + A1 + A2 + A3

    def mod_J_for_monomial(self, f, T, nth_roots_for_Rq, inv_nth_root_of_unity):
        """

        Input: - f an element in A can be written as f0 + f1*i + f2*j +f3*k represented with respect to monomial basis
               - T: the indices of T in E
               - n: the parameter that defines R
               - ninv: the inverse of n
               - nth_roots_for_Rq: pairs wich represent E
               - inv_nth_root_of_unity: the inv of the nth root of unity

        Output: f mod J reprsented with respect to monomial basis
        """
        n = self.smalln
        f0 = self.get_lagrange_coeff(f[0:n ** 2], nth_roots_for_Rq)
        f0 = self.mod_Q(f0, T)
        #     print("f0 after: ", f0)
        f0 = self.get_standard_coeff(f0, inv_nth_root_of_unity)

        f1 = self.get_lagrange_coeff(f[n ** 2:2 * n ** 2], nth_roots_for_Rq)
        #     print("f1: ", f1)
        f1 = self.mod_Q(f1, T)
        #     print("f1 after: ", f1)
        f1 = self.get_standard_coeff(f1, inv_nth_root_of_unity)

        f2 = self.get_lagrange_coeff(f[2 * n ** 2: 3 * n ** 2], nth_roots_for_Rq)
        f2 = self.mod_Q(f2, T)
        f2 = self.get_standard_coeff(f2, inv_nth_root_of_unity)

        f3 = self.get_lagrange_coeff(f[3 * n ** 2:4 * n ** 2],  nth_roots_for_Rq)
        f3 = self.mod_Q(f3, T)
        f3 = self.get_standard_coeff(f3, inv_nth_root_of_unity)
        F = f0 + f1 + f2 + f3

        return F

    def get_inverse_in_Aq(self, F, nth_roots_for_Rq, inv_nth_root_of_unity):
        """
        Input:  - F: the element to invert
                - qprime, n: parameters that defines the cryptosystem
                - ninv: n^{-1} in the field Fq
                - nth_roots_for_Rq, inv_nth_root_of_unity: defines the roots of unity and their inverses
        Output: - F^{-1}: the inverse of F in the quaternion A_q

        """
        # F =[0, 0, 0, -1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0]
        # print("norm: ", self.get_norm_bqtru(F))
        lagrange_coef_NF = self.get_lagrange_coeff(self.get_norm_bqtru(F), nth_roots_for_Rq)  ## N(F)
        # print("lagrange_coef", lagrange_coef_NF)
        NFinvF = self.get_coeff_wise_inv(lagrange_coef_NF)  ## (N(F))^{-1}
        # print("NFinv: ", NFinvF)
        standard_coeff_NFinvF = self.get_standard_coeff(NFinvF,  inv_nth_root_of_unity)  ## writting N(F) with respect to monomial bases.
        # print("standard", standard_coeff_NFinvF)
        fconjeguate = self.get_conjugate(F)  ## F
        # print("fconjeguate: ", fconjeguate)
        Finv = self.Rq_element_by_A(standard_coeff_NFinvF, fconjeguate)  ## F^{-1} = (N(F))^{-1} * F
        # print("F: ", F)
        # print("Finv: ", Finv)
        # print("get inverse in Aq: ", self.multiplication_in_A(F, Finv, self.q))
        return Finv

    def sample_W(self):
        """
        Input:  - qprime that defines the cryptosystem
        Output: - W = w0 + w1*i + w2*j + w3*k where w in L_q = (1,1/{Z_q})
        """
        while (True):
            w0 = self.FFq(randint(1, self.q - 1))
            w1 = self.FFq(randint(1, self.q - 1))
            w2 = self.FFq(randint(1, self.q - 1))
            w3 = self.FFq(randint(1, self.q - 1))

            norm = w0 ** 2 - w1 ** 2 - w2 ** 2 + w3 ** 2
            if norm != 0:
                break
        return [w0, w1, w2, w3]

    def w_sigma_monomial(self,w, sigma_monomial):
        """
        Input: - w = w0 + w1*i + w2*j + w3*k
               - sigma_monomial
        Output: -w*sigma
        """

        v0 = self.scalar_by_list(w[0], sigma_monomial)
        v1 = self.scalar_by_list(w[1], sigma_monomial)
        v2 = self.scalar_by_list(w[2], sigma_monomial)
        v3 = self.scalar_by_list(w[3], sigma_monomial)

        return v0 + v1 + v2 + v3

    def w_sigma_lagrange(self, w, sigma_lagrange):
        """
        Input: - w = w0 + w1*i + w2*j + w3*k
               - sigma_lagrange
        Output: - w*sigma represented in lagrange basis
        """
        return self.scalar_by_list(w[0], sigma_lagrange) + self.scalar_by_list(w[1], sigma_lagrange) + self.scalar_by_list(w[2],sigma_lagrange) + self.scalar_by_list(
            w[3], sigma_lagrange)


    """
    Input:  - F: an element of BQTRU ring has the form f0 + f1*i +f2*j +f3*k
    Output: - Rotations of F: all the quaternion rotations only [This function doesn't return the multivariate rotations] 
    """
    def get_quaternion_rotations(self, F):

        N = self.N
        f0 = F[0:N]
        f1 = F[N:2*N]
        f2 = F[2*N: 3*N]
        f3 = F[3*N:4*N]
        rotations_list = []

        rotations_list.append((0,F)) ## append F itself
        rotations_list.append((1, f1+f0+f3+f2))
        rotations_list.append((2, f2+self.scalar_by_list(-1, f3)+f0+self.scalar_by_list(-1, f1)))
        rotations_list.append((3, self.scalar_by_list(-1,f3)+f2+self.scalar_by_list(-1,f1)+f0))

        return rotations_list ##The rotations_list will contain all the rotations sotred as (index, element)
        
    """
    Input:   - F: the element to be rotated
             - the index of the rotation
    Output:  - the rotation to be returned according to the index
    The function supports the rotations due to the quaternion structure only for now.
    """
    def quaternion_rotation_by_index(self, F, index):
        N = self.N
        f0 = F[0:N]
        f1 = F[N:2 * N]
        f2 = F[2 * N: 3 * N]
        f3 = F[3 * N:4 * N]

        if index==0:
            return F
        elif index ==1:
            return f1+f0+f3+f2
        elif index==2:
            return f2+self.scalar_by_list(-1, f3)+f0+self.scalar_by_list(-1, f1)
        elif index==3:
            return  self.scalar_by_list(-1, f3)+f2+self.scalar_by_list(-1, f1)+f0
        else:
            raise ValueError("Not supported index!")


    """
    Input: an element of the underlying group-ring.
    The function builds the matrix of the group ring element and 
    return True if it's invertible, otherwise, it returns False.
    """
    def is_invertible_R_p(self,element):

        NF = self.get_norm_bqtru(element)
        mat = self.get_bivariate_mat(NF, IntegerModRing(self.p))
        if mat.is_invertible():
            return True
        return False

        # Fp_mat = self.element_to_matrix(element, self.FFp)
        # if Fp_mat.is_invertible():
        #     return True
        # return False


    """
        Input: seed
        Output: a key (f,g,h) where h = g*f^-1 mod(q, X^n-1)
    """
    def get_key(self, seed):
        if self.h!=None:  ##h has been initialized in the constructor
            return (None, None,self.h)
        print("seed: ", seed)

        bitseed = self.randomBitArray(self.sample_key_bits, self.seed_f)
        # print("tuple seed", tuple(bitseed))
        seedT = tuple(bitseed)
        if seedT in self.cache:
            return self.cache[seedT]

        else:
            print("weak instance:", self.weak_instance)
            f, g, Fp, Fq = self.sample_fg_bqtru(seed[0], seed[1], self.empty_fset, self.weak_instance) ####
            h = self.get_h_for_bqtru(Fq, g)

        self.cache[tuple(self.seed)] = (f, g, h)
        return (f, g, h)

    """
        Input: - decrypt=1 indicates we are calling the function during the decryption phase
        In the decryption phase, the person who decrypts know in advance the set s(he doesn't guess)
              
        Ouput: - matrix D that generates the Ideal Q = < sigma, q >
        
        This function is called only after calling the function guess_s so that the self.guessed_s is initialized
    """
    def build_lattice_D(self, decrypt=1):

        ns = self.N
        # V = RR ** ns
        vector_list = []
        if decrypt == 1:  ##no guessing, I am decrypting
            guessed_T_vectors  = []
            for i in self.s:
                guessed_T_vectors.append(self.lagrange_basis[i])
        else:
            guessed_T_vectors = self.get_geussed_T_vectors()

        for i in guessed_T_vectors:
            vector_list.append(vector(ZZ, i))  ###adding lagrange basis of T
            ##print(V.linear_dependence(vector_list, zeros='left'))

        for i in range(ns):
            qi = [0] * ns
            qi[i] = self.q
            vector_list.append(vector(ZZ, qi))

        mat = matrix(ZZ, vector_list)
        mat = mat.echelon_form()

        #     for i in removed_indices:
        #         qi = [0]*ns
        #         qi[i] = qprime
        #         vector_list.append(vector(ZZ, qi))
        #         temp = V.linear_dependence(vector_list, zeros='left')
        #         if (len(temp)==1):
        #             print("An error happend!!!")
        #             t = vector_list.pop()
        # #             print("on vector found not independent")
        # #             print("The vector: ", t)
        # #             print("The position: ", i)
        # #             print("The linear combination: ", temp)
        self.D = mat[0:ns, :]
        return mat[0:ns, :]

    def get_Mprivate(self, D):
        """
        Input:  - D: the basis of the lattice that defines the private ideal Q = <sigma, q> in R
        Output: - Aprivate: basis for the lattice the defines the private ideal J = <sigma, q> in A i.e., Q + Qi+ Qj + Qk
        The lattice has the form [D 0  0 0//  0  D 0 0 // 0 0 D 0// 0 0 0 D]
        """
        ns = D.nrows()  ## The dimension of D
        FF = IntegerRing()

        MS2 = MatrixSpace(FF, ns, ns)
        D = MS2.matrix(D)
        Zero = MS2.matrix(np.zeros([ns, ns]))

        # print(upper_right)

        return block_matrix(4, 4, [D, Zero, Zero, Zero, Zero, D, Zero, Zero, Zero, Zero, D, Zero, Zero, Zero, Zero, D])
    def get_h_for_bqtru(self, Fq, g):
        """
        This function calculates the public key for bqtru
        Input:  -Fq: the inverse of f
                -g: an element of the quaternion
        Output: -The public key h calculated as: h = Fq*g (mod q)
        """

        W = self.sample_W()
        V_monomial= self.w_sigma_monomial(W, self.sigma_monomial)
        # print("v monomial original: ", V_monomial)
        h = self.multiplication_in_A(Fq, g, self.q)
        # print("hstar original", h)
        h = self.addCN(h, V_monomial)  ##coefficients-wise addition
        # print("h orginal: ", h)
        self.V_monomial = V_monomial   ### We save the calculated value of v
        self.V_lagrange = self.get_lagrange_coeff_for_A(self.V_monomial, self.nth_roots_for_Rq )## We save the calculated value of v lagrange
        self.h = h
        return h

    """
        Input: seed
        Output: Coppersmith-Shamir basis for BQTRU Euclidean  lattice
        For option = 0, we return the lattice with no dimension reduction
        For option = 1, we return the lattice with one layer of dimension reduction
        Other options may be added later
        
        The message recovery attack is always considered to be a weak instance case.
    """
    def get_lattice(self, seed, update=False):

        if self.option == 0 :
            print("inside get lattice: option 0")
            if self.weak_instance:
                return self.get_lattice_option0(seed, update, self.attack_type)
            return  self.get_lattice_option0_not_weak(seed, update) ## get the lattice of the non-weak instance
        elif self.option == 1:
            print("inside get lattice: option 1")
            if self.weak_instance:
                return self.get_lattice_option1(seed, update, self.attack_type)

            return self.get_lattice_option1_not_weak(seed, update)
        else:
            raise ValueError("Not supported option for dimension reduction yet!")

    def get_lattice_option0(self, seed, update=False, attack_type=0):
        """
        Input:     - seed:    a certain seed for which we generate the lattice
                   - update:  False for the first time we build the lattice, True for updated values
                   - attack_type: 0 for key recovery attack, 1 for message recovery attack
        Output: Coppersmith-Shamir basis that define  for the Euclidean part of BQTRU
        with no dimension reduction.
        This function serves as a handler that decides whether to build the matrix for key recovery
        or for message recovery
        """
        if attack_type==0:
            return self.get_lattice_key_option0(seed,update)
        else:
            return  self.get_lattice_message_option0(seed, update) ##build the lattice for the message recovery



    def get_lattice_option1(self, seed, update=False, attack_type=0):
        """
         Input:     - seed:    a certain seed for which we generate the lattice
                    - update:  False for the first time we build the lattice, True for updated values
                    - attack_type: 0 for key recovery attack, 1 for message recovery attack
         Output: Coppersmith-Shamir basis that define  for the Euclidean part of BQTRU
         with one layer of  dimension reduction.
         This function serves as a handler that decides whether to build the matrix for key recovery
         or for message recovery
         """
        if attack_type == 0:
            return self.get_lattice_key_option1(seed, update)
        else:
            return self.get_lattice_message_option1(seed, update)  ##build the lattice for the message recovery

    def get_lattice_key_option0(self, seed, update=False):
        """
        Generate Coppersmith-Shamir lattice
        Input: - seed:    a certain seed for which we generate the lattice
               - update:  False for the first time we build the lattice, True for updated values

        Output: Coppersmith-Shamir lattice for the key recovery attack for the weak instance case
        """
        n = self.n
        q = self.q

        if update == False:  #first time to call the function
            f, g, h = self.get_key(seed)
            # we will save the values of f, g, h to access them from the attack file
            self.f = f
            self.g = g
            self.h = h



        self.guess_s() ##At this point, we will set the value of s
        if self.guessed_s == 'over':
            return 'failure'  ###failure means we have guessed all the values
        guessd_V = self.guess_V(self.guessed_s)
        self.hstar = self.get_hstar(guessd_V)
        h= self.hstar


        B = [[0] * 2 * n for i in range(2 * n)]
        for i in range(n):
            B[i][i] = q
        for i in range(n):
            B[n + i][n + i] = 1

        element_mat = self.element_to_matrix(h,self.FFq)  ##if attack_type =1, then tilde is 1, therefore
                                                                       ## we build the corresponding matrix
        for i in range(n):
            for j in range(n):
                B[n + i][j] = int(element_mat[i][j])

        return B


    def get_lattice_message_option0(self, seed, update=False):


        """
        This function performs embedding for the CVP to be converted into SVP

        The lattice to be built has the form

        [qI_n         & 0_n   &  0
        p(H^{tilde}_n) & I_n   &  0
        c             &  0    &  1
        ]
        Generate Coppersmith-Shamir lattice
        Input: - seed:    a certain seed for which we generate the lattice
               - update:  False for the first time we build the lattice, True for updated values

        Output: Coppersmith-Shamir lattice for the message recovery attack for the weak instance case
        According to my understanding this matrix should be [q & 0 // pH^{tilde} & I]
        where H^{tilde} as defined in the paper
        """
        n = self.n
        q = self.q
        # print("inside the not weak lattice")
        if update == False:
            f, g, h = self.get_key(seed)
            # we will save the values of f, g, h to access them from the attack file
            self.f = f
            self.g = g
            self.h = h  ## In the case of the message recover attack, the attack is built based on h and not h*

            #### get random message and encrypt the message

            message = self.sample_message()
            self.message = message ##save the message for checking purposes
            # print("The original message: ", message)
            c = self.encrypt(message, h)
            print("encrypt a message: ")
            self.ciphertext = c
            print(self.ciphertext)

            # print(c)
            # print("decrypting the message: \n")
            #
            # lattice = self.build_lattice_D(decrypt=1)
            # # print("The lattice is: ", lattice)
            # message_prime = self.decryption(c, lattice)
            # print("message prime: ", message_prime)
            # print("original message: ", message)
            #
            # if (message == message_prime):
            #     print("message decrypted successfully!! \n")

        h = self.h
        c = self.ciphertext
        B = [[0] * (2 * n + 1) for i in range(2 * n + 1)]
        for i in range(n):
            B[i][i] = q
        for i in range(n):
            B[n + i][n + i] = 1

        element_mat = self.element_to_matrix(h, self.FFq, tilde=1)  ##if attack_type =1, then tilde is 1, therefore
        ## we build the corresponding matrix
        for i in range(n):
            for j in range(n):
                B[n + i][j] = self.p*int(element_mat[i][j])
        # print(matrix(B).det())
        for i in range(n):
            B[2*n][i] = c[i]  ###fill the last line for the embedding

        B[2*n][2*n] = 1

        # print("The det: ")

        # print("compare to:  ", self.q**self.n)

        # print("B")
        # print(B)
        # print("matrix B")
        # print(matrix(B))
        # print(matrix(B).det())
        return B


    def get_lattice_option0_not_weak(self, seed, update= False):
        """
        Generate Coppersmith-Shamir lattice for the case of BQTRU not weak instance
        """
        n = self.n
        q = self.q
        print("inside the not weak lattice")
        if update == False:
            f, g, h = self.get_key(seed)
            # we will save the values of f, g, h to access them from the attack file
            self.f = f
            self.g = g
            self.h = h  ## In the non weak instance, we use h and not h*

        self.guess_s() ## At this point, we will set the value of s
        if self.guessed_s == 'over':
            return 'failure'  ##failure means we have guessed all the values

        ##For correctly guessed s we can build Mprivate = [D 0 0 0 // 0 D 0 0 // 0 0 D 0 // 0 0 0 D]
        Mprivate= self.get_Mprivate(self.build_lattice_D())

        B = [[0] * 2 * n for i in range(2 * n)]
        for i in range(n):
            for j in range(n):
                B[i][j] = Mprivate[i][j]
        for i in range(n):
            B[n + i][n + i] = 1

        element_mat = self.element_to_matrix(h,self.FFq)
        for i in range(n):
            for j in range(n):
                B[n + i][j] = int(element_mat[i][j])

        return B

    """
        Input: seed
        Output: Coppersmith-Shamir basis that define  for the Euclidean part of BQTRU 
        with one-layer dimension reduction.
        
        The lattice has the form 
        [qI 0  0   0 //  0 qI  0  0 //  H0+H1  H2+H3 I 0 // H2-H3   H0-H1 0 I]
    """

    def get_lattice_key_option1(self, seed, update = False, attack_type=0):

        n = self.n
        q = self.q

        if update == False:
            f, g, h = self.get_key(seed)
            # we will save the values of f, g, h to access them from the attack file to write the data in a file
            self.f = f
            self.g = g
            self.h = h


        self.guess_s()  ##At this point, we will set the value of s
        if self.guessed_s == 'over':
            return 'failure'  #### means I have guessed every set
        guessd_V = self.guess_V(self.guessed_s)

        ## get hstar to build the Eucledian lattice
        self.hstar = self.get_hstar(guessd_V)
        h = self.hstar
        B = [[0] * n for i in range(n)]
        half_n = int(n/2)
        for i in range(half_n):
            B[i][i] = q
        for i in range(half_n):
            B[half_n + i][half_n + i] = 1
        N = self.N  ###1/4 *n
        H0 = self.get_bivariate_mat(h[:N], self.FFq)
        H1 = self.get_bivariate_mat(h[N:2 * N], self.FFq)
        H2 = self.get_bivariate_mat(h[2 * N:3 * N],self.FFq)
        H3 = self.get_bivariate_mat(h[3 * N:4 * N], self.FFq)

        H0pH1 = H0 + H1
        H2pH3 = H2 + H3
        H2mH3 = H2 - H3
        H0mH1 = H0 - H1


        for i in range(N):
            for j in range(N):
                B[half_n+i][j] = int(H0pH1[i][j])   ##filling H0+H1

        for i in range(N):
            for j in range(N):
                B[half_n+i][N+j] = int(H2pH3[i][j])  ###filling H2+H3

        for i in range(N):
            for j in range(N):
                B[3*N+i][j] = int(H2mH3[i][j])     ### filling H2-H3

        for i in range(N):
            for j in range(N):
                B[3*N+i][N+j] = int(H0mH1[i][j])   ###filling H0-H1

        return B

    def get_lattice_message_option1(self, seed, update = False, attack_type=0):

        """

        Input: seed
        Output: Coppersmith-Shamir basis that defines the Euclidean part of BQTRU
        with one-layer dimension reduction.

        The function returns two lattices generated from the following basis:
        1.
        [qI 0_N  0_N   0_N  0// 0_N qI_N  0_N  0_N 0 0//  p(H0+H1)  p(H2-H3) I_N 0_N 0// p(H2+H3)   p(H0-H1) 0_N I_N 0// c0+c1 c2-c3 0_N 1_N 1]

        2.
        [qI 0_N  0_N   0_N // 0_N qI_N  0_N  0_N //  p(H0+H1)  p(H2-H3) I 0 // p(H2+H3)   p(H0-H1) 0 I// c2+c3 c0-c1 0_N 1]
        """

        n = self.n
        q = self.q
        # print("inside the not weak lattice")
        if update == False:
            f, g, h = self.get_key(seed)
            # we will save the values of f, g, h to access them from the attack file
            self.f = f
            self.g = g
            self.h = h  ## In the case of the message recover attack, the attack is built based on h and not h*

            #### get random message and encrypt the message

            message = self.sample_message()
            self.message = message  ##save the message for checking purposes
            # print("The original message: ", message)
            c = self.encrypt(message, h)
            # print("encrypt a message: ")
            self.ciphertext = c
            # print(self.ciphertext)
            ### The ciphertext is being interpreted as: (c0, c1, c2, c3)
        h = self.h
        c = self.ciphertext
        N = self.N  ###1/4 *n

        c0 = c[:N]
        c1 = c[N: 2*N]
        c2 = c[2*N:3*N]
        c3 = c[3*N:4*N]

        B = [[0] * n for i in range(n)]
        half_n = int(n / 2)
        for i in range(half_n):
            B[i][i] = q
        for i in range(half_n):
            B[half_n + i][half_n + i] = 1

        H0 = self.get_bivariate_mat(h[:N], self.FFq)
        H1 = self.get_bivariate_mat(h[N:2 * N], self.FFq)
        H2 = self.get_bivariate_mat(h[2 * N:3 * N], self.FFq)
        H3 = self.get_bivariate_mat(h[3 * N:4 * N], self.FFq)

        H0pH1 = H0 + H1
        H2pH3 = H2 + H3
        H2mH3 = H2 - H3
        H0mH1 = H0 - H1

        for i in range(N):
            for j in range(N):
                B[half_n + i][j] = self.p*int(H0pH1[i][j])  ##filling H0+H1

        for i in range(N):
            for j in range(N):
                B[half_n + i][N + j] = self.p*int(H2mH3[i][j])  ###filling H2-H3

        for i in range(N):
            for j in range(N):
                B[3 * N + i][j] = self.p*int(H2pH3[i][j])  ### filling H2+H3

        for i in range(N):
            for j in range(N):
                B[3 * N + i][N + j] = self.p*int(H0mH1[i][j])  ###filling H0-H1

        c0pc1 = self.addCN(c0,c1)         ### c0+c1
        c2mc3 = self.substractCN(c2,c3)   ### c2-c3
        c2pc3 = self.addCN(c2,c3)         ### c2+c3
        c0mc1 = self.substractCN(c0,c1)   ### c0-c1

        B1 = [[0] * (n+1) for i in range(n+1)]       ##B1
        B2 = [[0] * (n+1) for i in range(n+1)]       ##B2

        for i in range(n):
            for j in range(n):
                B1[i][j] = B[i][j]
                B2[i][j] = B[i][j] ##copy B1 to both of B1 and B2

        for i in range(N):
            B1[n][i] = c0pc1[i]
            B2[n][i] = c2pc3[i]

        for i in range(N):
            B1[n][N+i] = c2mc3[i]
            B2[n][N+i] = c0mc1[i]

        B1[n][n] = 1
        B2[n][n] = 1

        return (B1, B2)




    """
            Input: seed
            Output: Coppersmith-Shamir basis that defines the Euclidean part of BQTRU 
            with one-layer dimension reduction.

            The lattice has the form 
            [D 0  0   0 //  0 D  0  0 //  H0+H1  H2+H3 I 0 // H2-H3   H0-H1 0 I]
        """

    def get_lattice_option1_not_weak(self, seed, update=False):


        n = self.n
        q = self.q
        if update == False:

            f, g, h = self.get_key(seed)
            # we will save the values of f, g, h to access them from the attack file to write the data in a file
            self.f = f
            self.g = g
            self.h = h

        self.guess_s()  ##At this point, we will set the value of s
        D = self.build_lattice_D()  ##build D based on the guessed D

        # print("The det:", float(math.log(D.det(), self.q)))
        N = self.N  ###1/4 *n

        B = [[0] * n for i in range(n)]
        half_n = int(n / 2)
        for i in range(N):
            for j in range(N):
                B[i][j] = D[i][j]
                B[i + N][j + N] = D[i][j]

        for i in range(half_n):
            B[half_n + i][half_n + i] = 1

        H0 = self.get_bivariate_mat(h[:N], self.FFq)
        H1 = self.get_bivariate_mat(h[N:2 * N], self.FFq)
        H2 = self.get_bivariate_mat(h[2 * N:3 * N], self.FFq)
        H3 = self.get_bivariate_mat(h[3 * N:4 * N], self.FFq)

        H0pH1 = H0 + H1
        H2pH3 = H2 + H3
        H2mH3 = H2 - H3
        H0mH1 = H0 - H1

        for i in range(N):
            for j in range(N):
                B[half_n + i][j] = int(H0pH1[i][j])  ##filling H0+H1

        for i in range(N):
            for j in range(N):
                B[half_n + i][N + j] = int(H2pH3[i][j])  ###filling H2+H3

        for i in range(N):
            for j in range(N):
                B[3 * N + i][j] = int(H2mH3[i][j])  ### filling H2-H3

        for i in range(N):
            for j in range(N):
                B[3 * N + i][N + j] = int(H0mH1[i][j])  ###filling H0-H1

        return B



    """
    Sample a message
    """

    def sample_message(self):

        random = Random()
        if self.seedm == None:
            self.seedm = randint(0, 2**64)

        random.seed(self.seedm)

        ############ uncomment this block for random ternary messages ##############
        # seed_for_m = [random.randrange(2) for i in range(self.sample_iid)]
        # seqlen = int(len(seed_for_m)/4)
        # m0 = self.ternary(seed_for_m[0:seqlen])
        # m1 = self.ternary(seed_for_m[seqlen: 2*seqlen])
        # m2 = self.ternary(seed_for_m[2*seqlen: 3*seqlen])
        # m3 = self.ternary(seed_for_m[3*seqlen: 4*seqlen])
        # self.message = m0+m1+m2+m3
        # return m0+m1+m2+m3

        ################## uncomment this block for messages coming from the space L(dm, dm)
        seed_for_m = [random.randrange(2) for i in range(self.sample_fixed_type)]
        seqlen = int(len(seed_for_m) / 4)
        m0 = self.fixed_type(seed_for_m[0:seqlen], self.d, self.d)
        m1 = self.fixed_type(seed_for_m[seqlen: 2 * seqlen], self.d, self.d)
        m2 = self.fixed_type(seed_for_m[2 * seqlen: 3 * seqlen], self.d, self.d)
        m3 = self.fixed_type(seed_for_m[3 * seqlen: 4 * seqlen], self.d, self.d)
        self.message = m0+m1+m2+m3

        return m0+m1+m2+m3


    """
      Input: a message to encrypt
             h: the public key
      Output: the encrypted message
      """

    def encrypt(self, message, h):

        random = Random()
        if self.seedr==None:
            self.seedr = randint(0, 2 ** 64)

        random.seed(self.seedr)
        # if self.seedr!=None:
        #     random.seed(self.seedr)   #### generate r according to the already determined seedr
        # else:
        #     random.seed(randint(0, 2 ** 64))

        seed_for_r = [random.randrange(2) for i in range(self.sample_fixed_type)]
        seqlen = int(len(seed_for_r)/4)
        r0 = self.fixed_type(seed_for_r[0:seqlen], self.d, self.d)
        r1 = self.fixed_type(seed_for_r[seqlen: 2*seqlen], self.d, self.d)
        r2 = self.fixed_type(seed_for_r[2*seqlen: 3*seqlen], self.d, self.d)
        r3 = self.fixed_type(seed_for_r[3*seqlen: 4*seqlen], self.d, self.d)
        r = r0+r1+r2+r3
        # print("r: ", r)
        self.r = r ##save it so you can access the value of r in the attack file and write it there
        e1 = self.multiplication_in_A(h, r, self.q)
        prh = list(np.multiply(self.p, e1))
        e = list(np.add(prh, message))
        self.ciphertext = e
        return e

    def reduce_lattice(self, B, option=0, beta=None, float_type="mpfr"):
        """
        Input:  - B: a basis for the lattice
                - option: to apply for lattice reduction
                      0: no lattice reduction (default)
                      1: LLL
                      2: BKZ with appropriate blocksize beta
               - float_type: the default is d for larger lattices mpfr can be used for precision
        """

        B = IntegerMatrix.from_matrix(B)
        if option == 0:
            return B
        else:
            ## Create GSO object, I am creating a copy of B and pass it to the function so that I can do the check for M.U*B = M.B
            M = GSO.Mat(copy(B), float_type=float_type, U=IntegerMatrix.identity(B.nrows, int_type=B.int_type),
                        UinvT=IntegerMatrix.identity(B.nrows, int_type=B.int_type))
            ## no need to update now
            bkz = BKZReduction(M)  ###create an object of BKZreduction

            if option == 1:
                bkz.lll_obj()

            else:
                ## apply the reduction algorithm with the required blocksize
                par = BKZ_FPYLLL.Param(beta, strategies=BKZ_FPYLLL.DEFAULT_STRATEGY, max_loops=8)  ##Parameters
                bkz(par)

            return M.B
    ########## decryption ################

    def get_CVP_in_D(self, D, target):
        """
        Input:  - D: the private lattice where CVP to be solved, we have Mprivate = [D 0 0 0 // 0 D 0 0 // 0 0 D 0//0 0 0 D]
                - target  : the target vector
        Output: - a vector in the lattice that is the closed to the target vector
        """
        reduced_D = self.reduce_lattice(D, 1)  ## reduced lattice (you can apply different options)
        n = self.smalln ## the value of n
        closet_vector1 = CVP.closest_vector(reduced_D, target[0:n**2])
        closet_vector2 = CVP.closest_vector(reduced_D, target[n**2:2*n**2])
        closet_vector3 = CVP.closest_vector(reduced_D, target[2*n**2:3*n**2])
        closet_vector4 = CVP.closest_vector(reduced_D, target[3*n**2:4*n**2])

        return closet_vector1 + closet_vector2 + closet_vector3 + closet_vector4

    def decryption(self, c,  lattice):
        """
        Input: -  c:         the encrypted message
               -  lattice:   the lattice D
        """


        FC = self.multiplication_in_A(self.f, c, self.q)
        target = vector(ZZ, self.center_lift_form(FC, self.q))
        B = self.get_CVP_in_D(lattice, target)  ##Find the closed vector to FC in Aprivate where Aprivate is the lattice and FC is the target vector
        # print("B: ", B)
        V = self.addCN(self.center_lift_form(FC, self.q),
                          self.scalar_by_list(-1, B))  ### FC-B: a quaternion with integer coefficients
        V = self.center_lift_form(V, self.q)
        # print("V: ", V)
        # print("norm: ", float(get_norm(V)))
        #     print("Fp: ", Fp)
        #     print("V: ", V)
        mprime = self.multiplication_in_A(self.Fp, V, self.p)
        mprime = self.center_lift_form(mprime, self.p)
        return mprime

    """
    The function returns the value of f
    """
    def get_f(self):
        return self.f

    """
    The function returns the value of g
    """
    def get_g(self):
        return  self.g


    """
    The function returns the public key h
    """
    def get_h(self):
        return self.h


    """
    The function returns the message
    """
    def get_message(self):
        return  self.message

    """
    The function returns the ciphertext
    """
    def get_ciphertext(self):
        return  self.ciphertext



    """
    The function returns the saved r
    """
    def get_r(self):

        return self.r

    def guess_s(self):

        """
        The functions guesses the positions of the T points that defines sigma
        """
        if self.guess == False:
            self.guessed_s = self.s  ###for now, no guessing
                                     ## We are just returning the guessed value
        else:
            self.guessed_s = self.get_guess_as_set() #if we need to guess, we go one by one
            self.guess_next()

        # print("guessed s: ", self.guessed_s)


    def guess_V(self, guessed_s):
        """
        Input: -  guessed_s: the position of T that defines sigma
        Output: - V_monomial: if the guessed positions are correct

        """
        print("guessed_s", guessed_s)
        # print("self.N", self.N)
        # print("self.g inside guess_V", self.g)
        # eval_G = []
        # for i in range(self.N):
        #     if i in guessed_s:
        #         point = self.nth_roots_for_Rq[i]
        #         eval_G.append(self.evaluate_two_variables(self.g[0:self.N], point))
        #
        # for i in range(self.N):
        #     if i in guessed_s:
        #         point = self.nth_roots_for_Rq[i]
        #         eval_G.append(self.evaluate_two_variables(self.g[self.N:2*self.N], point))
        #
        # for i in range(self.N):
        #     if i in guessed_s:
        #         point = self.nth_roots_for_Rq[i]
        #         eval_G.append(self.evaluate_two_variables(self.g[2*self.N:3*self.N], point))
        #
        #
        # for i in range(self.N):
        #     if i in guessed_s:
        #         point = self.nth_roots_for_Rq[i]
        #         eval_G.append(self.evaluate_two_variables(self.g[3*self.N:4*self.N], point))
        #
        #
        # print("eval G: ", eval_G)

        self.guessed_s = guessed_s ##save the guessed s

        V_guessed_lagrange = []
        for i in range(self.N):
            if i in guessed_s:
                point = self.nth_roots_for_Rq[i]
                V_guessed_lagrange.append(self.evaluate_two_variables(self.h[0:self.N],point))
            else:
                V_guessed_lagrange.append(0)

        for i in range(self.N):
            if i in guessed_s:
                point = self.nth_roots_for_Rq[i]
                V_guessed_lagrange.append(self.evaluate_two_variables(self.h[self.N:2*self.N], point))
            else:
                V_guessed_lagrange.append(0)

        for i in range(self.N):
            if i in guessed_s:
                point = self.nth_roots_for_Rq[i]
                V_guessed_lagrange.append(self.evaluate_two_variables(self.h[2*self.N:3*self.N], point))
            else:
                V_guessed_lagrange.append(0)

        for i in range(self.N):
            if i in guessed_s:
                point = self.nth_roots_for_Rq[i]
                V_guessed_lagrange.append(self.evaluate_two_variables(self.h[3*self.N:4*self.N], point))
            else:
                V_guessed_lagrange.append(0)

        # print("V_guessed lagrange: ", V_guessed_lagrange)
        V_guessed_monomial0 = self.get_standard_coeff(V_guessed_lagrange[0:self.N], self.inv_nth_roots_of_unity)
        V_guessed_monomial1 = self.get_standard_coeff(V_guessed_lagrange[self.N:2*self.N], self.inv_nth_roots_of_unity)
        V_guessed_monomial2 = self.get_standard_coeff(V_guessed_lagrange[2*self.N:3*self.N], self.inv_nth_roots_of_unity)
        V_guessed_monomial3 = self.get_standard_coeff(V_guessed_lagrange[3*self.N:4*self.N], self.inv_nth_roots_of_unity)
        V_guessed_monomial = V_guessed_monomial0+V_guessed_monomial1+V_guessed_monomial2+V_guessed_monomial3
        print("V guessed monomial: ", V_guessed_monomial)
        print("V lagrange again: ", self.get_lagrange_coeff_for_A(V_guessed_monomial, self.nth_roots_for_Rq))
        return V_guessed_monomial


    def get_hstar(self, guessed_V_monomial):
        """
        This function guess v and then return hstar = h -v

        To guess v, you need to guess T different positions and from these positions, you get v
        For trivial F, we have F^{-1}*F = 1 mod q
        We know that H = F^{-1}*G+v mod q ==> by guessing T positions, we can get
        H(a1, b1), H(a2, b2),......H(aT, bT) that gives v , because the evaluation at these points on G is zero

        For now, we will get the saved value of v, later on we can search for it
        """
        # print("h down: ", self.h)
        hstar = self.Modulo(self.substractCN(self.h, guessed_V_monomial), self.q)

        # print("hstar down: ", hstar)
        # print("checking the calculations are correct: it should come to be the ternary g")
        # print("f down: ", self.f)
        # print("fq down: ", self.Fq)
        # print("check f*fq: ", self.multiplication_in_A(self.f, self.Fq, self.q))
        t= self.multiplication_in_A(self.f, hstar, self.q)
        # print("The multiplication result: ", self.center_lift_form(t))
        # print("g: ", self.g)
        # assert(self.center_lift_form(t)==self.g)
        return hstar

    """
    Input: seed
    Output: Reduction lattices
            - in the case of cyclic group, it's one lattice
            - For dihedral group, two lattices of dimension 2*n
    """

    def get_reduction_lattices(self, seed):
        n = self.n
        q = self.q
        f, g, h = self.get_key(seed)




    """
    This function to be deleted 
    We are using this function just for testing purposes.
    """
    def checking_invertibility(self, t, f_set_empty=True):
        """

        :return:
        """

        self.update_seed(True)
        f_seed = self.seed_f
        for i in range(t):

            # print("s for G: ", s)
            seqlen = int(self.sample_fixed_type / 4)
            while (True):
                f_bits = self.randomBitArray(self.sample_fixed_type, f_seed)

                ### sample F
                f0 = self.fixed_type(f_bits[0:seqlen], self.d + 1, self.d)
                f1 = self.fixed_type(f_bits[seqlen:2 * seqlen], self.d, self.d)
                f2 = self.fixed_type(f_bits[2 * seqlen:3 * seqlen], self.d, self.d)
                f3 = self.fixed_type(f_bits[3 * seqlen:4 * seqlen], self.d, self.d)
                f = f0 + f1 + f2 + f3
                NF = self.get_norm_bqtru(f)
                fset = self.get_T(NF)
                if f_set_empty:
                    condition = len(fset) == 0  ##sampling f with empty set
                else:
                    condition = (len(fset) >= 1)
                if (condition):  ### and len(fset)>=1
                    print("The condition of the set cardinality is satisfied !")
                    mat = self.get_bivariate_mat(NF, IntegerModRing(self.p))
                    if mat.is_invertible():
                        Fp = mat.inverse()[0]
                        print("The invertibility condition is satisfied through method 1 !")
                        Fp = self.Rq_element_by_A(Fp, self.get_conjugate(f))
                        print("The inverse through the first approach", Fp)


                        #############
                        mat = self.element_to_matrix(f, IntegerModRing(self.p))
                        if  self.is_invertible_R_p(f):
                            print("The element is inveritble through the second approach")
                            Fp2 = mat.inverse()[0]
                            print("The inverse through the second approach", Fp2)

                        break
                self.update_seed(True)  ##update seed_f
                print("updating seed f")
                f_seed = self.seed_f


def main():

    # n = 7*4
    # q = 63
    # print(is_it_ternary([1, -1, 0, -1, 1, 0, 0, -1, 0, 0, 1, 1, 0, 0]))
    # for i in range(100):
    #
    #     seed = randint(0,2**64)
    #     keygen = NTRUKeyGenerator(n, q,seed=seed, qtru= True)
    #     keyseed = keygen.newSeed()
    #
    #     # f,g,h= keygen.get_key(keyseed)
    #     # print(h)
    #     # print(keygen.multiply_for_QTRU(f,h,q))
    #     # print(g)
    #     print(keygen.get_lattice(keyseed))


    n = 4*5**2
    q = 241
    counter = 0
    for i in range(1):
        seed_f = randint(0, 2 ** 64)


        seed_g0 = randint(0, 2 ** 64)
        seed_g1 = randint(0, 2 ** 64)
        seed_g2 = randint(0, 2 ** 64)
        seed_g3 = randint(0, 2 ** 64)

        seed_g = (seed_g0, seed_g1, seed_g2, seed_g3)
        seed = (seed_f, seed_g)

        #### empty_fset = True, we will test our program for this option
        #### There are some instances where the decryption does not work successfully, and these instances happen
        #### when the subset size is greater than 6
        keygen = NTRUKeyGenerator(n, q, seed=seed, option=0, empty_fset= True, weak_instance=True)
        g,f,h = keygen.get_key(seed)
        print("The public key: ", h)
        #print(keygen.element_to_matrix(h, FF=IntegerRing(), tilde=1))
        message = keygen.sample_message()
        print("The orginal message: ", message)
        c= keygen.encrypt(message,h)
        print("encrypt a message: ")
        print(c)
        print("decrypting the message: \n")

        lattice = keygen.build_lattice_D(decrypt=1)
        # print("The lattice is: ", lattice)
        message_prime = keygen.decryption(c,lattice)
        print("message prime: ", message_prime)
        print("original message: ", message)

        if(message==message_prime):
            print("message decrypted successfully!! \n")
            counter = counter+1
        else:
            print("message did not decrypt successfully for the set {}".format(keygen.s))

    L = keygen.get_lattice_message_option0(seed, update=True)


    print("counter is: ", counter)





        # t = keygen.get_guess_as_set()
        # while (True):
        #     print(t)
        #     print("-----------------------------------------")
        #     keygen.guess_next()
        #     t = keygen.get_guess_as_set()
        #     if t=='over':
        #         break
        #
        # print(t)
        # print("-----------------------------------------")
        # print("counting: ", keygen.counting)
        # print(keygen.get_bivariate_mat([1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15, 16]))

        # (f,g, Fp, Fq) = keygen.sample_fg_bqtru(seed_f, seed_g)
        # print(f)
        # #print(keygen.get_BQTRU_mat([1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15, 16]))
        #
        # # keygen.checking_invertibility(10,True)
        #
        # keygen.get_lattice((seed_f, seed_g))
        # # print(keygen.get_seed())
        #
        # # # print(keygen.V_monomial)
        # # # print(keygen.V_lagrange)
        # print("guessing")
        # guessed_V = keygen.guess_V(keygen.s)
        # print(guessed_V==keygen.V_monomial)
        # hstar = keygen.get_hstar(guessed_V)
        # print("hstar: ", hstar)
        # f = keygen.get_f()
        # g = keygen.get_g()
        # print("f: ", f)
        # print("g: ", g)
        # checked_g =  keygen.center_lift_form(keygen.multiplication_in_A(f,hstar, keygen.q))
        # assert (checked_g== g)
        # one = keygen.center_lift_form(keygen.multiplication_in_A(f, keygen.Fq, keygen.q))
        # print("one: ", one)
        # g_supposed = keygen.center_lift_form(keygen.multiplication_in_A(f, keygen.hstar, keygen.q))
        # print("g_supposed: ", g_supposed)
        # assert (g == g_supposed)
        # print(keygen.multiplication_in_A(one, keygen.h, keygen.q))
        # print(keygen.h)


    # n = 11
    # q = 64
    #
    # seed = randint(0, 2**64)
    # keygen = NTRUKeyGenerator(n, q, seed=seed, option=0)
    # lattice = keygen.get_lattice(seed)
    # print(lattice)
    # print(get_dual(lattice))

if __name__ == "__main__":
    main()