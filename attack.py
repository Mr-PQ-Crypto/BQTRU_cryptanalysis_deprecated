import time
from fpylll import BKZ as BKZ_FPYLLL, GSO, IntegerMatrix, FPLLL
from fpylll.tools.quality import basis_quality
from fpylll.algorithms.bkz2 import BKZReduction
from random import randint
import json

from httplib2 import re_unsafe

import keygen
from keygen import NTRUKeyGenerator
from utils import (get_key_norm, get_norm, is_it_zero, is_it_ternary, is_it_pm_2, add_vectors_with_centerlifting,
                   substract_vectors_with_centerlifting, divide_by_2, run_all, parse_args, dump_blocksize_for_group,
                   dump_seed, dump_blocksize_for_message_attack, sliceBasis, scalar_multiply, addCN, substractCN, is_it_boundedby)

FPLLL.set_precision(120)
"""
Implementation of the lattice reduction for both cyclic group(NTRU) and dihedral group(DiTRU)


How to:
   

Examples:

for dihderal group: 
    python attack.py 14  -q=128 --verbose=True --dump=True --group="dihedral" --h="[115, 42, 117, 108, 73, 3, 53, 29, 108, 34, 72, 5, 36, 101]"
    python attack.py 89 --verbose=True --group="cyclic" --dump=True --h="[403, 317, 342, 288, 441, 171, 102, 79, 33, 92, 295, 146, 160, 287, 480, 264, 167, 164, 396, 216, 493, 219, 351, 163, 140, 384, 62, 290, 218, 410, 151, 101, 369, 8, 398, 118, 231, 152, 428, 233, 172, 55, 307, 480, 58, 86, 469, 45, 10, 317, 121, 399, 234, 26, 108, 498, 325, 234, 37, 228, 456, 6, 371, 475, 310, 3, 22, 378, 419, 196, 93, 158, 124, 409, 286, 187, 216, 490, 71, 69, 122, 261, 388, 405, 71, 397, 343, 470, 337]" --filename="test" --bkz_betas="3:15"
"""

class Attack:
    """
    Attack on NTRU over cyclic group and dihedral group.
    We find the first block size that find the key for equivalent instances
    It includes the block size needed to find a {ternary key, non_ternary}

    Examples:

    for dihedral group:
            python attack.py 14  -q=128 --verbose=True --group="dihedral" --h="[115, 42, 117, 108, 73, 3, 53, 29, 108, 34, 72, 5, 36, 101]"
    """

    def __init__(self, params):

        self.n = params['n']  # The order of the group
        self.q = params['q']  # The used modulo
        self.option = params['option'] # option = 0 default means no dimension reduction
                                       # option = 1 one layer of dimension reduction
                                       # option = 2,... to be discussed later
        self.group = 'bqtru'
        self.file_tag = params['file_tag'] # File tag
        self.empty_fset = params['empty_fset']
        self.weak_instance = params['weak_instance']

        self.seed_f  = params['seed_f']
        self.seed_g0 = params['seed_g0']
        self.seed_g1 = params['seed_g1']
        self.seed_g2= params['seed_g2']
        self.seed_g3 = params['seed_g3']
        self.seed = (self.seed_f, (self.seed_g0, self.seed_g1, self.seed_g2, self.seed_g3) )


        self.seed_r  = params['seed_r']
        self.seed_m  = params['seed_m']

        self.guess = params['guess'] ## if True means do guessing otherwise assume you have guessed corrctly and do
                                    # lattice reduction directly
        #self.nsamples = params['nsamples']  # number of NTRU samples
        self.blocksizes = params['blocksizes']  # range [a,b] with starting blocksize a, last blocksize b-1
        self.ntours = params['tours']  # number of bkz tours
        self.nthreads = params['threads']  # number of threads
        self.verbose = params['verbose']  # verbose mode
        self.filename = params['filename'] ##file name
        self.dump     = params['dump']
        self.message    = None
        self.ciphertext = None
        self.attack_type = params['attack_type'] #attack type, 0 for key recovery attack and 1 for memory recovery attack
        keyGenSuccess = False
        if params['h'] != None:
            self.h = json.loads(params['h'])  # for user-input h, we don't know f,g.
            self.f = None
            self.g = None
        else:
            self.h = None

        self.generator = NTRUKeyGenerator(self.n, self.q, self.seed, self.h, self.option, self.empty_fset, self.weak_instance, self.guess, self.attack_type)
        if self.attack_type ==1:  ###for message attack
            self.generator.set_seed_r(self.seed_r)
            self.generator.set_seed_m(self.seed_m)

        # while not keyGenSuccess:
        #     try:
        #         self.keyseed = self.generator.newSeed()
        #         #print("seed: ", self.seed)
        #         #print(self.generator.get_key(self.keyseed))
        #         self.lattice = self.generator.get_lattice(self.keyseed)
        #         # self.lattice contains a tuple: as (lattice, plus_lattice, minus_lattice)
        #         # in the case of the cyclic group both of plus_lattice and minus_lattice are None
        #
        #         keyGenSuccess = True
        #     except:
        #         self.seed = self.seed+1
        #         self.generator.update_seed(self.seed)

        self.lattice = self.generator.get_lattice(self.seed)  ##return the reduced dimension lattice if option=1

        self.seed    = self.generator.get_seed()
        print("updated seed after  calling self.lattice", self.seed)
        ### At this point we have the lattice of the element either in Z_qC_n or Z_qD_n
        ## Also the object generator is saving f,g corresponding to the lattice.
        ## We are getting these values in order to save them in the file
        self.f = self.generator.get_f()
        self.g = self.generator.get_g()
        self.h = self.generator.get_h()


        if self.attack_type == 1:##the attack is message recovery
            self.message    =      self.generator.get_message()
            self.ciphertext =  self.generator.get_ciphertext()
            self.r          = self.generator.get_r()
            self.seed_r     = self.generator.get_seedr()
            self.seed_m     = self.generator.get_seedm()
            # print("The ciphertext value", self.ciphertext)

        if self.attack_type==1 and self.option==1:
            self.dim = len(self.lattice[0])
        else:
            self.dim = len(self.lattice)

        # print("dim: ", self.dim)
        # print("one vector: ", len(self.lattice[0]))
        # self.dim = self.n * 2  # For cyclic group the lattice dimension for the reduction is two times the order

        if self.option ==0:
            factor = 4
        else:
            factor =2

        self.threshold = factor * get_key_norm(self.n)

        if self.dim <= 178:
            self.float_type = "long double"
        else:
            self.float_type = "mpfr"

        if self.attack_type==1 and self.option==1:
            # The case of message attack and dimension reduction


            self.basis1 = IntegerMatrix.from_matrix(self.lattice[0], int_type="long")
            self.M1 = GSO.Mat(self.basis1, float_type=self.float_type,
                             U=IntegerMatrix.identity(self.basis1.nrows, int_type=self.basis1.int_type),
                             UinvT=IntegerMatrix.identity(self.basis1.nrows, int_type=self.basis1.int_type))

            self.basis2 = IntegerMatrix.from_matrix(self.lattice[1], int_type="long")
            self.M2 = GSO.Mat(self.basis2, float_type=self.float_type,
                              U=IntegerMatrix.identity(self.basis2.nrows, int_type=self.basis2.int_type),
                              UinvT=IntegerMatrix.identity(self.basis2.nrows, int_type=self.basis2.int_type))

        else:
            self.basis = IntegerMatrix.from_matrix(self.lattice, int_type="long")
            self.M = GSO.Mat(self.basis, float_type=self.float_type,
                             U=IntegerMatrix.identity(self.basis.nrows, int_type=self.basis.int_type),
                             UinvT=IntegerMatrix.identity(self.basis.nrows, int_type=self.basis.int_type))





    def __call__(self):
        # print("f, which we are looking for: ", self.f)
        # g = self.generator.multiplication_in_A(self.f, self.generator.hstar, self.q)
        # print("g: should be: ", self.generator.center_lift_form(g))
        ###dump the seed before
        if self.dump:
            if self.attack_type ==0:
                dump_seed(self.seed, self.group, self.filename, self.file_tag)
            else:
                dump_seed(self.seed, self.group, self.filename, self.file_tag, self.attack_type, (self.seed_m, self.seed_r))
        self.progressive_search()  # call the function that retrieves the key


    def call_update(self, previous_time=0.0):

        """
        When we are guessing the values, in case we are not getting the solution
          - update the set of guessing, then build the lattice using the same options of the original call
          - call the progressive search again in the range and see if we are getting the solution
        """
        self.lattice = self.generator.get_lattice(self.seed, True)

        self.basis = IntegerMatrix.from_matrix(self.lattice, int_type="long")
        self.M = GSO.Mat(self.basis, float_type=self.float_type,
                         U=IntegerMatrix.identity(self.basis.nrows, int_type=self.basis.int_type),
                         UinvT=IntegerMatrix.identity(self.basis.nrows, int_type=self.basis.int_type))

        self.progressive_search(previous_time)



    def check_for_layer0(self, found, target):

        """
        Upon a reduced basis of a lattice for BQTRU
        This function call either check_for_layer0 for the message
        or check_for_layer0 for the key
        Input: - found: True if what are we looking for got found, False otherwise
               - target: passes the found keys/messages
        """

        if self.attack_type == 0:
            out = self.check_for_layer0_key(found, target)
        else:
            out = self.check_for_layer0_message(found, target)

        return  out

    def check_for_layer0_message(self, message_found, message):

        """
        Upon a reduced basis of a lattice for BQTRU
        check if a vector that look like the encoded message has been found
        message_found: refers if the message has been found or not
        message: returns the message
        """
        ### message coming as [[],[]]
        norms = {}
        B = self.M.B
        for i in range(self.dim):
            norms[i] = B[i].norm()
        sorted_norms = sorted(norms.items(), key=lambda x: x[1])
        for i in range(self.dim):
            if sorted_norms[i][1]>self.threshold:
                if not message_found: ##threshold is the ||(m, r)||
                  return  "failure"  ##if we are meeting vectors with norm greater than the target (m, r)
                else:
                    return message

            mr = list(B[sorted_norms[i][0]])
            # print("checking the vector mr", mr)
            if not(is_it_zero(mr[:self.n])) and is_it_ternary(mr):
                message_found = True
                if mr[self.dim-1]==-1:
                    mr = scalar_multiply(-1, mr)   ###The tuple has the form (m, -r, 1)
                print("possible message found")
                message[0].append(mr[0:self.n])  ### getting m (we append all possible messages)
                message[1].append(scalar_multiply(-1,mr[self.n:]))  ### getting r

        return  "failure"



    def check_for_layer0_key(self, keys_found_tuple, keys):
        """
        Upon a reduced basis of a lattice for BQTRU
        check for the ternary/non-ternary key.
        keys_found_tuple: a tuple refers if the (non-ternary-found, ternary-found).
        keys: it is passed as (None, None) initially, then will save the results of the found keys over rounds
        Output: a tuple (k1, k2) where ki itself is a tuple as ([f,g], norm).
        if no key is returned, returns "failure".
        """

        key1_found = keys_found_tuple[0]
        key2_found = keys_found_tuple[1]
        k1 = keys[0]  ##(key, norm)
        k2 = keys[1]  ##(key, norm)
        norms = {}
        # print(self.M.B)
        B = self.M.B
        for i in range(self.dim):
            norms[i] = B[i].norm()
        sorted_norms = sorted(norms.items(), key=lambda x: x[1])
        for i in range(self.dim):
            if sorted_norms[i][1] > self.threshold:
                if key1_found or key2_found:
                    return (k1, k2)
                return "failure"

            fg = list(B[sorted_norms[i][0]])
            f = fg[self.n:]
            g = fg[:self.n]
            if not is_it_zero(g):
                # print("checking vector : \n")
                # print(f+g)
                # print("check it gives small value: ",
                #       self.generator.multiplication_in_A(f, self.generator.hstar, self.q))
                if not key1_found and self.generator.is_invertible_R_p(f):
                    k1 = (fg, sorted_norms[i][1])  # (key, its norm)

                    #print("check it gives small value: ", self.generator.multiplication_in_A(f, self.generator.hstar, self.q))
                    key1_found = True

            if not key2_found and is_it_ternary(fg):
                if self.generator.is_invertible_R_p(f):
                    k2 = (fg, sorted_norms[i][1])  # (key, its norm)
                    #print("check it gives small value: ", self.generator.multiplication_in_A(f, self.generator.hstar, self.q))
                    key2_found = True

            if key1_found and key2_found:
                return (k1, k2)

        return "failure"

    def check_for_layer1(self, found, target):

        """
        Upon a reduced basis of a lattice for BQTRU
        This function call either check_for_layer1 for the message
        or check_for_layer1 for the key
        Input: - found: True if what are we looking for got found, False otherwise
               - target: passes the found keys/messages
        """

        if self.attack_type == 0:
            out = self.check_for_layer1_key(found, target)
        else:
            out = self.check_for_layer1_message(found, target)

        return out

    def check_for_layer1_key(self, keys_found_tuple, keys):
        """
                Upon a reduced basis of a lattice for BQTRU,
                check for the ternary/non-ternary key.

                keys_found_tuple: a tuple refers of the (non_ternary found, ternary found).
                keys: it is passed as (None, None) initially, then will save the results of the found keys over rounds.
                Output: a tuple (k1, k2) where ki itself is a tuple as ([f,g], norm).
                if no key is found, returns "failure".
        """
        key1_found = keys_found_tuple[0]
        key2_found = keys_found_tuple[1]
        k1 = keys[0]
        k2 = keys[1]
        N = int(self.dim / 4)  # every vector is divided into four halves
        norms= {}
        B = self.M.B ## Will be the reduced-dimensional lattice
        for i in range(self.dim):
            norms[i] = B[i].norm()
        sorted_norms = sorted(norms.items(), key=lambda x: x[1])

        for i in range(self.dim):
            if sorted_norms[i][1] > self.threshold:
                if key1_found or key2_found:
                    return (k1, k2)
                return "failure"
            t = list(B[sorted_norms[i][0]])
            #print("N: ", N)
            g0 = t[0: N]
            g1 = t[N: 2*N]
            f0 = t[2*N: 3*N]
            f1 = t[3*N: 4*N]

            if (not is_it_zero(g0+g1)) and (not is_it_zero(f0+f1)):
                for j in range(self.dim):
                    if sorted_norms[j][1] > self.threshold:
                        break
                    t = list(B[sorted_norms[j][0]])
                    g0p = t[0:N]
                    g1p = t[N:2*N]
                    f0p = t[2*N: 3*N]
                    f1p = t[3*N: 4*N]

                    if (not is_it_zero(g0p+g1p)) and (not is_it_zero(f0p+f1p)):

                        fp0 = add_vectors_with_centerlifting(f0, f1p, N, self.q)
                        fp1 = substract_vectors_with_centerlifting(f0, f1p, N, self.q)
                        fp3 = add_vectors_with_centerlifting(f1, f0p, N, self.q)
                        fp4 = substract_vectors_with_centerlifting(f1, f0p, N, self.q)

                        gp0 = add_vectors_with_centerlifting(g0, g1p, N, self.q)
                        gp1 = substract_vectors_with_centerlifting(g0, g1p, N, self.q)
                        gp3 = add_vectors_with_centerlifting(g1, g0p, N, self.q)
                        gp4 = substract_vectors_with_centerlifting(g1, g0p, N, self.q)
                        F = fp0 + fp1 + fp3 + fp4 # concatenating (fp0, fp1)
                        G = gp0 + gp1  + gp3 + gp4 # concatenating  (gp0, gp1)


                        Flist = self.generator.get_quaternion_rotations(F)
                        for Fi in Flist:
                            index = Fi[0]
                            F     = Fi[1]
                            G = self.generator.quaternion_rotation_by_index(F, index)
                            if not key1_found and self.generator.is_invertible_R_p(F):
                                # print("invertible")
                                #print("(f0,g0): ", f0+g0)
                                #print("first vecor norm: ", get_norm(g0 + f0))
                                #print("(f1,g1): ", f1+g1)
                                #print("second vecor norm: ", get_norm(g1 + f1))
                                k1 = (F + G, get_norm(F + G))
                                key1_found = True

                            if not key2_found and is_it_pm_2(F + G):
                                F = divide_by_2(F)
                                G = divide_by_2(G)

                                if self.generator.is_invertible_R_p(F):
                                    k2 = (F + G, get_norm((F + G)))
                                    key2_found = True

                                if key1_found and key2_found:
                                    return (k1, k2)
        #print("reached here")
        return "failure"


    def check_for_layer1_message(self, message_found, message):
       """
        Upon a reduced basis of BQTRU,
        check for the message that can be lifted from the defined homomorphisms
        Input: message_found: True or False indicates whether the message found or not
               message: a list defined as [[],[]], a list of two lists, where the first list is used to store the
               possible found messages and the second list is used to save possible found -r

        Output: a list where ki itself is a tuple as ([f,g], norm).
        if no key is found, returns "failure".
       """

       # print("inside the check function! ")
       N = int(self.dim/4)
       norms1 = {}
       norms2 = {}
       B1 = self.M1.B  ## we have two different lattices in the case of one layer of the message attack
       B2 = self.M2.B

       # print("best inside the function: ", B1[0])
       # print("best inside the functipm: ", B2[0])
       for i in range(self.dim):
           norms1[i] = B1[i].norm()
           norms2[i] = B2[i].norm()
       sorted_norms1 = sorted(norms1.items(), key=lambda x: x[1])
       sorted_norms2 = sorted(norms2.items(), key=lambda x: x[1])

       for i in range(self.dim):
           if sorted_norms1[i][1] > self.threshold:
               if message_found:
                   return message
               return "failure"
           t = list(B1[sorted_norms1[i][0]])
           if t[self.dim-1] == -1:
               t = scalar_multiply(-1,t) ### we multiply by -1 to get [(m_0+m_1), (m_2-m_3), -(r_0+r_1), -(r_2-r_3), 1]
           m0pm1  = t[0:N]
           m2mm3  = t[N:2*N]
           mr0pr1 = t[2*N:3*N]
           mr2mr3 = t[3*N:4*N]

           if ((not is_it_zero(m0pm1 + m2mm3)) and (not is_it_zero(mr0pr1 + mr2mr3)) and
                   is_it_boundedby(m0pm1 + m2mm3,2) and is_it_boundedby(mr0pr1 + mr2mr3,2)):

               for j in range(self.dim):
                   if sorted_norms2[j][1] > self.threshold:
                       break
                   t = list(B2[sorted_norms2[j][0]])
                   if t[self.dim-1] == -1:
                       t = scalar_multiply(-1,t) ## we multiply by -1 to get ##[(m_2+m_3), (m_0-m_1), -(r_2+r_3), -(r_0-r_1), 1]
                   m2pm3  = t[0:N]
                   m0mm1  = t[N:2 * N]
                   mr2pr3 = t[2 * N: 3 * N]
                   mr0mr1 = t[3 * N: 4 * N]

                   if ((not is_it_zero(m2pm3 + m0mm1)) and (not is_it_zero(mr2mr3 + mr0mr1)) and
                       is_it_boundedby(m2pm3 + m0mm1, 2) and is_it_boundedby(mr2mr3 + mr0mr1,2)):
                       m02  = addCN(m0pm1, m0mm1) ##getting 2*m0
                       m12  = substractCN(m0pm1, m0mm1)## getting 2*m1
                       m22  = addCN(m2mm3, m2pm3) ##getting 2*m2
                       m32  = substractCN(m2pm3, m2mm3) ###getting 2*m3


                       mr02 = addCN(mr0pr1, mr0mr1) ### getting -2*r0
                       r12  = substractCN(mr0pr1, mr0mr1) ###getting 2*r1
                       mr22  = addCN(mr2mr3, mr2pr3) ###getting -2*r2
                       r32  = substractCN(mr2pr3, mr2mr3) ###getting 2*r3


                       if (is_it_pm_2(m02) and is_it_pm_2(m12) and is_it_pm_2(m12)
                           and is_it_pm_2(m22) and is_it_pm_2(m32) and is_it_pm_2(mr02)
                       and is_it_pm_2(r12) and is_it_pm_2(mr22) and is_it_pm_2(r32)):
                           m0 = divide_by_2(m02)
                           m1 = divide_by_2(m12)
                           m2 = divide_by_2(m22)
                           m3 = divide_by_2(m32)

                           r0 = divide_by_2(scalar_multiply(-1, mr02))
                           r1 = divide_by_2(r12)
                           r2 = divide_by_2(scalar_multiply(-1,mr22))
                           r3 = divide_by_2(r32)

                           m = m0+m1+m2+m3
                           r = r0+r1+r2+r3
                           message_found = True
                           message[0].append(m)
                           message[1].append(r)

       return "failure"









    def progressive_search(self, previous_time=0.0):
        """
         Apply reduction algorithm with increased block sizes and return the minimum blocksize
         that retrieves both a non-ternary and ternary keys
         previous time = 0 for no previous call
                       otherwise, we add the time recursively
        """


        print("inside progressive search: ")

        continue_bkz = False
        stop_bkz = False
        key1_found = False  # The non ternary key found?.
        key2_found = False # The ternary key found?.
        key_found = (key1_found, key2_found)
        key1 = (None, None) ## (The non ternary key, its norm)
        key2 = (None, None) ## (The ternary key, its norm)
        message_list = [[],[]]
        key_tuple = [key1, key2]
        message_found = False
        if self.attack_type==0:
            target_tuple = key_tuple
            target_found = key_found
        else:
            target_tuple = message_list
            target_found = message_found
        beta = [0]*2 ## block size needed to retrieve the (non-ternary key, ternary key)

        T0_global = time.time()

        # print("before: ", self.M.B)

        if self.attack_type==1 and self.option==1:
            ##In the case of the message attack with one layer of reduction
            self.bkz1 = BKZReduction(self.M1)
            self.bkz1.lll_obj()

            self.bkz2 = BKZReduction(self.M2)
            self.bkz2.lll_obj()
        else:
            self.bkz = BKZReduction(self.M)
            self.bkz.lll_obj()
        if self.option==0:
            target_tuple = self.check_for_layer0(target_found, target_tuple)
        else:
            target_tuple = self.check_for_layer1(target_found, target_tuple)

        if self.verbose:
        # print("group:", self.group)
            fmt = "{'initial LLL applied: underlying group':'%8s', 'total walltime': %.3f}"
            print(fmt % (self.group, time.time() - T0_global))
        if target_tuple == "failure":
            continue_bkz = True
            if self.attack_type ==0:
                target_tuple = [key1, key2]
            else:
                target_tuple = [[],[]]
            if self.verbose:
                print("failure")
        else:
            if self.attack_type == 0:
                if (target_tuple[0][0]!=None):
                    if self.verbose:
                        print("(Non ternary key, its norm)", target_tuple[0])
                    key1_found = True
                    beta[0] = 2
                if target_tuple[1][0]!=None:
                    if self.verbose:
                        print(("(Ternary key ,its norm)", target_tuple[1]))
                        # print("check it's giving small g, for the ternary key: ", self.generator.multiplication_in_A(key_tuple[1], self.generator.hstar, self.q))
                    key2_found = True
                    beta[1] = 2

                if not (key1_found and key2_found):
                    print("continue_bkz")
                    continue_bkz = True
            else:  ##message attack
                 if (len(target_tuple[0])!=0):
                     if self.verbose:
                         beta[0] = 2
                         message_found = True
                         print("possible messages")
                         for mess in target_tuple[0]:
                             print(mess)
                 else:
                     continue_bkz = True

        if continue_bkz:
            for blocksize in self.blocksizes: ##apply bkz with increasing block size

                T0_local = time.time()
                if self.verbose:
                    print("New round with block size: ", blocksize)
                    if self.attack_type==1 and self.option==1:  ##message attack with one-layer of reduction
                        quarter = self.generator.N
                        m0 = self.message[0:quarter]
                        m1 = self.message[quarter:2*quarter]
                        m2 = self.message[2*quarter: 3*quarter]
                        m3 = self.message[3*quarter: 4*quarter]

                        r0 = self.r[0:quarter]
                        r1 = self.r[quarter:2*quarter]
                        r2 = self.r[2*quarter: 3*quarter]
                        r3 = self.r[3*quarter: 4*quarter]

                        # m0pm1 = addCN(m0, m1)
                        # m2m3  = substractCN(m2, m3)
                        # r0pr1 = addCN(r0, r1)
                        # mr0pr1 = scalar_multiply(-1, r0pr1)
                        # r2mr3 = substractCN(r2, r3)
                        # mr2mr3 = scalar_multiply(-1, r2mr3)

                        # print("[ [(m_0+m_1), (m_2-m_3), -(r_0+r_1), -(r_2-r_3), 1] ]")
                        # print(m0pm1+m2m3+mr0pr1+mr2mr3+[1])

                        print("best vector in  phi+ \n")
                        print(self.M1.B[0])

                        m2pm3 = addCN(m2, m3)
                        m0mm1 = substractCN(m0, m1)
                        r2pr3 = addCN(r2, r3)
                        mr2pr3 = scalar_multiply(-1, r2pr3)
                        r0mr1 = substractCN(r0, r1)
                        mr0mr1 = scalar_multiply(-1, r0mr1)

                        # print("[(m_2 + m_3), (m_0 - m_1), -(r_2 + r_3), -(r_0 - r_1), 1]")
                        # print(m2pm3+m0mm1+mr2pr3+mr0mr1+[1])
                        print("best vector in phi- \n")
                        print(self.M2.B[0])

                        # print("M1 B")
                        # print(self.M1.B)
                        # print("M2 B")
                        # print(self.M2.B)
                    else:
                        print("best vector in the reduced lattice \n")
                        print(self.M.B[0])


                for t in range(self.ntours):  # runs BKZ tours
                     par = BKZ_FPYLLL.Param(blocksize,
                                            strategies=BKZ_FPYLLL.DEFAULT_STRATEGY,
                                            max_loops=8)

                     if self.attack_type==1 and self.option ==1:  ##for message attack with one layer of reduction
                         self.bkz1(par)
                         self.bkz2(par)
                     else:
                         self.bkz(par)  ##apply bkz with the required parameters the check again for key
                     #print("new tour")
                     if self.option==0:
                         target_tuple = self.check_for_layer0(target_found, target_tuple)
                     else:
                         target_tuple = self.check_for_layer1(target_found, target_tuple)


                     #print("key tuple after calling the function: ", key_tuple[0], "   ", key_tuple[1])
                     if target_tuple =="failure":
                         if self.attack_type == 0:
                             target_tuple = [key1, key2]
                         else:
                             target_tuple = [[],[]]
                         if self.verbose:
                             print("failure")

                     else:
                         if self.attack_type ==0:

                             if target_tuple[0][0] != None and not  key1_found:
                                 if self.verbose:
                                     print("(Non ternary key, its norm)", target_tuple[0])
                                 key1_found = True
                                 beta[0] = blocksize
                             if target_tuple[1][0] != None:
                                 if self.verbose:
                                     print(("(Ternary key ,its norm)", target_tuple[1]))
                                     print("key_tuple: ", target_tuple[1][0][self.n:2*self.n])
                                     #print("check it's giving small g, for the ternary key: ", self.generator.multiplication_in_A(key_tuple[1][0][self.n:2*self.n], self.generator.hstar, self.q))
                                 key2_found = True
                                 beta[1] = blocksize

                             if key1_found and key2_found:
                                 stop_bkz = True
                                 break

                         else: ### message attack
                             if (len(target_tuple[0]) != 0):
                                 message_found = True
                                 beta[0] = blocksize
                                 if self.verbose:
                                     print("possible messages")
                                     for mess in target_tuple[0]:
                                         print(mess)
                             if message_found:
                                 stop_bkz = True
                     if stop_bkz:
                         break
                     fmt = "{'BKZ 'beta': %2d, underlying group':'%8s', 'total walltime': %.3f}"
                     if self.verbose:
                        print(fmt % (blocksize,self.group, time.time() - T0_local))
                if stop_bkz:
                    break
        if self.attack_type == 0:
            if self.verbose:
                print("Block size needed to find (non-ternary key, ternary key) is ({},{})".format(beta[0], beta[1]))

            if beta[0] == 0 or beta[1]==0:  # I have reached to the end of the range and no solution came
                previous_time = previous_time+ (time.time()-T0_global)
                if self.verbose:
                    fmt = "{ 'Total walltime till this trial' %.3f}"
                    print(fmt % (previous_time))
                return  self.call_update(previous_time)

            if self.dump:
                #print("non ternary key: ", key_tuple[0])
                #print("ternary key: ", key_tuple[1])
                dump_seed(self.seed,self.group,self.filename, self.file_tag)
                dump_blocksize_for_group(self.f, self.g, self.h, target_tuple, beta, self.filename, self.file_tag, self.group, self.generator.guessed_s, self.generator.V_lagrange, self.generator.V_monomial, (time.time()-T0_global)+previous_time)

        else: ##message attack
            if self.verbose:
                print("Block size needed to find possible message is {}".format(beta[0]))

            if self.dump:
                # dump_seed(self.seed, self.group, self.filename, self.file_tag, self.attack_type, (self.seed_m, self.seed_r))
                dump_blocksize_for_message_attack(self.h,self.ciphertext, self.message, self.r,  target_tuple[0],target_tuple[1] ,  beta[0], self.filename, self.file_tag, self.group, (time.time()-T0_global))
                ###Tomorow here
def main_call(params):
    attack_inst = Attack(params)
    return attack_inst()


if __name__ == '__main__':

    all_parsed_params = parse_args()
    run_all(main_call, all_parsed_params)
