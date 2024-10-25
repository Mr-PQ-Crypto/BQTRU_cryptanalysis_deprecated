import argparse
import six
import copy
import sys, os
import re
from math import sqrt, floor, ceil, log, exp, log2
from datetime import datetime
import random
from random import randint
from multiprocessing import Pool
from sage.all import Matrix, ZZ
from random import randint
from fpylll.util import gaussian_heuristic
from fpylll import IntegerMatrix, LLL
from mpmath import mp, erfc as mperfc
from scipy.special import  erfcinv
import csv
from sympy import  nextprime
mp.dps = 1200
error_rate = 2**-64
import warnings



max_beta = 150 # maximal dimension we run SVP in
cross_over_beta = 70

def get_dual(basis):
    """
    Input: basis vector of a lattice
    Output: the dual matrix of the basis

    """
    B = Matrix(basis)
    B = B.pseudoinverse().transpose()
    B = B.denominator() * B
    B = B.change_ring(ZZ)

    basis = []
    for v in B:
        basis.append(v.list())

    return basis


def swapLeftRight(basis):
  """
    On input of a basis matrix B = (B_1||B_2),
    this function computes (B_2||B_1).
  """
  d = int( len(basis[0]) /2)

  for i in range(len(basis)):

    left = basis[i][:d]
    right = basis[i][d:]
    basis[i] = right + left

  return basis


def projection(basis, projectLeft):
  """
    On input of a basis matrix B = (B_1||B_2) and projectLeft = True,
    this function projects B_1 orthogonally against the all one vector.
    Then scales the matrix, such that it is integral.
    If projectLeft = False, then the projection is applied to B_2.
  """
  d = int( len(basis[0]) /2)

  if not projectLeft:
    basis = swapLeftRight(basis)

  for i, v in enumerate(basis):
    v_left = v[:d]
    v_right = v[d:]

    sum_left = sum(v_left)

    for j in range(d):
      v_left[j] = d*v_left[j] - sum_left
      v_right[j] = d*v_right[j]

    basis[i] = v_left + v_right

  if not projectLeft:
    basis = swapLeftRight(basis)

  return basis


def removeLinearDependencies(basis):
  """
    Removes linear dependencies using LLL:
  """

  t = []

  for b in basis:
      t.append(list(b))

  d = len( basis[0] )

  # print("d: ", d)
  B = IntegerMatrix(d, d)

  # B.set_matrix(Matrix(t))
  B = IntegerMatrix.from_matrix(basis)
  B = LLL.reduction(B)
  while list(B[0]) == B.ncols*[0]:
    B = B[1:]

  basis = [ [ 0 for _ in range(B.ncols) ] for _ in range(B.nrows) ]
  B.to_matrix( basis )

  return basis

def sliceBasis(basis,projectLeft=True):
  """
    On input of a basis matrix for
    the Coppersmith-Shamir lattice (lattype=classic) or
    the projected cylcotomic lattice (lattype=phi_projected),
    this function computes a basis for the lattices with additional hints by design,
    i.e., classic_slive or phi_projected_slice
    as introduced in Section 5.3. of our paper.
  """
  basis = get_dual(basis)
  basis = projection(basis,projectLeft)
  basis = removeLinearDependencies(basis)
  basis = get_dual(basis)
  print("slice basis is over")
  return basis

def remove_zeros(B):
    """
    removes zero rows from matrix B
    """
    cr = 0
    for i in range(B.nrows):
      if not B[i].is_zero(): cr+=1

    B_ = [0]*cr
    cr = 0
    for i in range(B.nrows):
      if not B[i].is_zero():
        B_[cr] = B[i]
        cr+=1

    return IntegerMatrix.from_matrix(B_, int_type="long")

def testAndMakeDir(path):
  if not os.path.isdir(path):
      os.makedirs(path)


def get_norm(vector):
    """
    Calculates and returns the norm of the vector.
    """
    # print("inside the function: \n")
    # print("iside vector: ", vector)

    s = 0
    for v in vector:
        s += v * v
    # print("inside the norm: ", np.sqrt(s))
    return sqrt(s)

def get_key_norm(n):
    """
    Input: n the order of the group
    Output: the norm of the key of the form (f,g)
    where f in T(d+1,d) and g in T(d,d).
    """

    d = int(n/3)
    return sqrt(4*d+1)

def is_it_ternary(l):
    """
    Input: a list l.
    Check if the list is ternary and returns True otherwise returns True
    """
    for i in l:
        if i!=1 and i!=-1 and i!=0 :
            return False
    return True

def is_it_boundedby(l, bound):
    """

    :param l: list
    :param bound: bound to check if all the elements in the list are bounded by
    :return: True if all the elements in the list are bounded by the bound and False otherwise
    """
    for i in l:
        if(abs(i)>bound):
            return False
        return True

def is_it_zero(l):
    """
    Input: list l
    Output: True if the list entries are all zeros, False otherwise
    """
    for i in l:
        if i!=0:
            return False
    return True

def get_q_no_error(d,p):
    """
    The function returns the value of q that gives no decryption failure for variant of NTRU
    that has: h = gf^-1
    Input: d = int(order of the group/3)
           p usually 3
    """
    value= p*(6*d+1)
    q= 2**(len(bin(value))-2)
    return q



def is_it_pm_2(l):
    """
    Input: a list.
    Return True if all entries are two, minus two, or zeros.
    Otherwise: False.
    """
    for i in l:
        if i!=2 and i!=-2 and i!=0 :
            return False
    return True


def divide_by_2(l):
    """
    Input: a list of {2,-2,0}
        divide the coefficients by 2 and return the resultant list.

    """
    for i in range(len(l)):
        if l[i]>0:
            l[i] =1
        elif l[i]<0:
            l[i] = -1
    return l
def add_vectors_with_centerlifting(l1, l2, n, q):
    """
    Coefficients-wise adding of coefficients of the correspondence vectors
    with centrelifting
    Input: l1, l2 two lists representing two vectors in the lattice.
    n: the vector length.
    Output: the resultant vector after adding them.
    """
    res = [0]*n
    for i in range(n):
        res[i] = (l1[i]+l2[i])%q
        if res[i]>int(q/2):
            res[i] = res[i]-q
    return res


def substract_vectors_with_centerlifting(l1, l2, n, q):
    """
    Coefficients-wise adding of coefficients of the correspondence vectors
    with centrelifting
    Input: l1, l2 two lists representing two vectors in the lattice.
        n: the vector length
    Output: the resultant vector after adding them.
    """
    res = [0] * n
    for i in range(n):
        res[i] = (l1[i] - l2[i]) % q
        if res[i] > int(q / 2):
            res[i] = res[i] - q
    return res


def addCNmodq(a, b, mod):
    """
    Input: a,b: two elements in C_N
           mod: the mod of adding
    Output: c: adding the two vectors
    """
    n = len(a)
    c = [0] * n
    for i in range(n):
        c[i] = (a[i] + b[i]) % mod
    return c

def addCN(a, b):
    """
    Input: a,b: two elements in C_N

    Output: c: adding the two vectors coefficients-wise without modulo
    """
    n = len(a)
    c = [0] * n
    for i in range(n):
        c[i] = (a[i] + b[i])
    return c


def substractCN( a, b):
    """
    Input: a,b two elements in C_N
    Output: adding the corresponding coefficients.
    """
    n = len(a)
    c = [0] * n
    for i in range(n):
        c[i] = a[i] - b[i]
    return c

def scalar_multiply(sc, li):
    """
    Input: sc: a scalar value to be multiplied by the list li
    Output: the list li multiplied by the scalar li
    """
    output_list = [sc * ele for ele in li]
    return output_list
def dump_seed(seed, group,filename, filetag, attack_type=0, seed_mr=None):
    """
    Input: the seed, the file name, and the filetag, attack_type
            attack_type =0 for key recovery attack and 1 for message recovery
    Output: write the seed to the file to later add the trails and the betas.
    """

    org_seed = seed
    if seed_mr!=None:
        org_seed = seed + seed_mr ###adding
    # seed = seed - (seed % 10 ** 8)
    if attack_type==0:
        path = "keys_dumps/" + group + "/" +"seeds/"
    else:
        path = "messages_dumps/" + group + "/" + "seeds/"
    testAndMakeDir(path)

    filename += "_" + str(filetag) + ".txt"
    with open(path + filename, "a+") as f:
        print("seed: ",org_seed, file=f)



def process_sample(sample):
    """

    Input: a sample as  a list of info
    a sample can be a list holding infor like: ['f', 'g', 'key norm', 'h', 'k1 (non-ternary)', 'k1-norm', 'k2 (ternary)', 'k2-norm', 'beta1', 'beta2', 'total time (seconds)']

    The function processes the sample and returns a list of strings values corresponding to the values in the list.

    """
    l = []

    for i in range(len(sample)):

        s = str(sample[i])
        string_to_write = ""
        for i in range(len(s)):
            if (s[i] != "]" and s[i] != "["):
                string_to_write += s[i]
        l.append(string_to_write)
    return l


def create_file(file_tage, group,filename,attack_type=0):
    """
    Input: the seed, the group and the file name, and attack_type(0 for key recovery attack and 1 for message recovey)
    The function creates a file with the specified path and write the header to the file.
    the header = [f, g, norm, h, f1_prime, f1_norm, f_prime2, f2_norm, beta1, beta2, total_time]
    The function writes the header into a csv file, if not already existed.
    """
    #print("file created: ")
    print("file tag: ", file_tage)
    # seed = seed - (seed % 10 ** 8)
    if attack_type ==0:
        path = "keys_dumps/" +group + "/records/"
        header = ['f', 'g', 'key norm', 'h', 'k1 (non-ternary)', 'k1-norm', 'k2 (ternary)', 'k2-norm', 'beta1', 'beta2',
                  'guessed T', 'guessed v(lagrange)', 'guessed v(monomial)', 'total time (seconds)']

    else:
        path = "messages_dumps/" +group + "/records/"
        header = ['h', 'c', 'original: m', 'original: r', 'retrieved: m', 'retrieved: r', 'beta', 'total time (seconds)']
    testAndMakeDir(path)


    filename += "_" + str(file_tage) + ".csv"
    isExisting = os.path.exists(path+filename)
    if not isExisting:
        with open(path + filename, "w", newline='') as wfl:
           csvwriter = csv.writer(wfl, delimiter=',')
           csvwriter.writerow([val for val in header])

def dump_blocksize_for_group(f,g,h,key_tuple, beta, filename,file_tag, group, guessed_s, guessed_v_lagrange, guessed_v_monomial, total_time):
  """
  Input:
         f,g,h: the original key (f,g): is the private key and h is the corresponding public key.
         key_tuple: it is an array key_tuple[0] = (non-ternary-key, its norm) and key_tuple[1] = (ternary key, its norm)
         key_tuple can be a string "failure" if we couldn't find the key.
         beta: an array: beta[0]: the blocksize needed to find the non-ternary key and beta[1]: the blocksize needed to
         find the ternary key.
         filename: the file name (recommended to be "n_q")
         group: bqtru
         guessed_s: the guessed value of s
         guessed_v_lagrange: the guessed value of v with respect to lagrange basis
         guessed_v_monomial: the guessed value of v with respect to monomial basis
         seed: the seed.
         total)time: the total time the attack took to find the keys
  """
  sample = []
  sample.append(f)
  sample.append(g)
  if f== None:
      sample.append("NA")
  else:
      sample.append(get_norm(f+g))
  sample.append(h)

  if key_tuple =="failure":
      sample.append("failure") ## not able to find the non-ternary key
      sample.append("NA")  ## no norm
      sample.append("failure") ## not able to find the ternary key
      sample.append("NA") ## no norm
      sample.append("failure for tried betas") ## beta for the non ternary key
      sample.append("failure for tried betas") ## beta for the ternary key
  else:
      sample.append(key_tuple[0][0]) ## non ternary key
      sample.append(key_tuple[0][1]) ## non-ternary norm
      sample.append(key_tuple[1][0]) ## ternary key
      sample.append(key_tuple[1][1]) ## ternary norm
      sample.append(beta[0]) ## beta1 needed to find the non-ternary key
      sample.append(beta[1]) ## beta2 needed to find the ternary key

  sample.append(guessed_s)
  sample.append(guessed_v_lagrange)
  sample.append(guessed_v_monomial)

  sample.append(total_time)

  # seed = seed - (seed % 10 ** 8)
  path = "keys_dumps/" + group + "/records/"
  testAndMakeDir(path)

  to_write = process_sample(sample)
  filename += "_" + str(file_tag) + ".csv"
  with open(path + filename, "a+", newline='') as wfl:
      csvwriter = csv.writer(wfl, delimiter=',')
      csvwriter.writerow([val for val in to_write])
      # print( str(beta1) + "\t" + str(beta2) + "\t" + str(total_time) + "\t"+ datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file = f )



def dump_blocksize_for_message_attack(h,c,original_message,original_r, m,r ,  beta, filename, file_tag, group, totaltime):
    """

    :param h: the public key
    :param c: the ciphertext
    :param original_message: the original message
    :param original_r: r sampled in the encryption process
    :param m: the retrieved message
    :param r: the retrieved r
    :param beta: the blocksize needed to retrieve the message
    :param filename: the file name
    :param file_tag: file tage
    :param group: bqtru
    :param totaltime: the totaltime to do the attack
    :return: create file and save the mentioned info in it
    """

    sample  = []
    sample.append(h)
    sample.append(c)
    sample.append(original_message)
    sample.append(original_r)
    if len(m)==0:
        sample.append("failure") ##for m
        sample.append("failure") ##for r
        sample.append("NA") # beta
    else:
        sample.append(m)
        sample.append(r)
        sample.append(beta) ##beta is passed as a single value for the message recovery attack
    sample.append(totaltime)
    path = "messages_dumps/" + group + "/records/"
    testAndMakeDir(path)

    to_write = process_sample(sample)
    filename += "_" + str(file_tag) + ".csv"
    with open(path + filename, "a+", newline='') as wfl:
        csvwriter = csv.writer(wfl, delimiter=',')
        csvwriter.writerow([val for val in to_write])

def rough_estimate_on_betas(n,q):
    """
    use the output of this function in case the use did not provide us with blocksizes
    """

    if n<150: beta_low = 10
    else: beta_low = floor( 0.28*4.*n/( (log(q)/log(n))**2 + 1))
    return list(range(beta_low, max_beta))




def parse_args():
    parser = argparse.ArgumentParser(description='Parse NTRU attack params.')

    #main parameters
    parser.add_argument('n', type=int, help="ring dimension")
    parser.add_argument('file_tag', type=str, help="the file tag") ##we make the file tag mandatory because we can't generate the file uniquely based on the seed due to the complicated seeds associated with bqtru
    parser.add_argument('-q', type=int, dest="q", default=None, help="NTRU modulus")
    parser.add_argument('--nsamples', type=int, default=None, dest="nsamples", help="Number of samples/rows of rot(h) used")

    parser.add_argument('--seed_f', type=int, dest="seed_f", default=None,   help="randomness seed for f (bqtru)")
    parser.add_argument('--seed_g0', type=int, dest="seed_g0", default=None, help="randomness seed for g0 (bqtru)")
    parser.add_argument('--seed_g1', type=int, dest="seed_g1", default=None, help="randomness seed for g1 (bqtru)")
    parser.add_argument('--seed_g2', type=int, dest="seed_g2", default=None, help="randomness seed for g2 (bqtru)")
    parser.add_argument('--seed_g3', type=int, dest="seed_g3", default=None, help="randomness seed for g3 (bqtru)")


    parser.add_argument('--seed_r', type=int, dest="seed_r", default=None, help="randomness seed for r")
    parser.add_argument('--seed_m', type=int, dest="seed_m", default=None, help="randomness seed fpr m")
    parser.add_argument('--option', type=int, dest="option", default=0, help="option=0 for no dimension reduction, option=1 for one layer of dimension reduction, ...etc")
    parser.add_argument('--h', dest="h", default=None, help="Uses given input as h, instead of creating a random instance.")
    parser.add_argument('--empty_fset', dest="empty_fset", default=True, help="if True, we generate a key with empty fset")
    parser.add_argument('--dump', dest='dump', default=False, help="flag to dump intermediate bases")
    # number of runs, number of threads
    parser.add_argument('-t', '--trials', type=int, dest="trials", default=1,
                        help="number of experiments to run per dimension")
    parser.add_argument('-w', '--workers', type=int, dest="workers", default=1,
                        help="number of parallel experiments to run")
    parser.add_argument('--threads', type=int, dest="threads", default = 1, help="number of threads used by 1 worker")
    parser.add_argument('--weak_instance', dest="weak_instance", default=True, help="if True, we generate the T set as mentioned in the paper T= {(a_i, b_i) in E s.t. g_i(a_i, b_i) = 0 } ")

    parser.add_argument('--bkz_betas', type=str, dest="blocksizes", default=None, help="bkz block sizes as string of the form: min_beta:max_beta:step")
    parser.add_argument('--bkz_tours', type=int, dest="tours", default=8, help="number of tours of bkz reduction")

    parser.add_argument('--guess', dest="guess", default=False, help="guess the positions of T!")
    parser.add_argument('--verbose', dest="verbose", default=False, help="verbosity")
    parser.add_argument('--dry-run', dest="dry_run", default=False,
                        help="Show parameters that would be used but don't run any actual experiments.")
    parser.add_argument('--show-defaults', dest="show_defaults", action='store_true',
                        help="Show default parameters and exit.")

    parser.add_argument('--attack_type', type=int, dest="attack_type", default=0, help="key recovery {0} or message recovery {1}")

    parser.add_argument('--filename', dest='filename', default=None, help="prefix of the dump filenames")



    args, unknown = parser.parse_known_args()


    fmt = "{key:%ds}: {value}"%20

    if len(unknown)>0:
        print('Parameters', unknown, 'are not recognized and will be ignored')

    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    if args.show_defaults:
        for k, v in six.iteritems(all_defaults):
            print(fmt.format(key=k, value=v))
        exit(0)

    all_params = check_parsed_params(vars(args))

    if args.dry_run:
        for k, v in six.iteritems(all_params):
            print(fmt.format(key=k, value=v))
        exit(0)


    return all_params

def check_parsed_params(params):

    if params['verbose']=="False" or params['verbose']==False:
        params['verbose'] = False
    else:
        params['verbose'] = True

    if params['weak_instance'] == "False" or params['weak_instance'] == False:
        params['weak_instance'] = False
    else:
        params['weak_instance'] = True

    if params['empty_fset'] == "False" or params['empty_fset']==False:
        params['empty_fset'] = False
    else:
        params['empty_fset'] = True


    if params['guess'] == "False" or params['guess'] ==False:

        params['guess'] = False
    else:
        params['guess'] = True


    if params['nsamples']==None: params['nsamples'] = params['n']
    else: assert(params['nsamples'] > 0 and params['nsamples']<=params['n'])

    if params['blocksizes'] ==None:
        params['blocksizes'] = rough_estimate_on_betas(params['n'], params['q'])
    else: params['blocksizes'] = eval("range(%s)" % re.sub(":", ",", params['blocksizes']))

    assert(len(params['blocksizes'])>0)

    if params['seed_f'] == None:
        params['seed_f'] = randint(0, 2 ** 64)

    if params['seed_g0'] == None:
        params['seed_g0'] = randint(0, 2 ** 64)
    if params['seed_g1'] == None:
        params['seed_g1'] = randint(0, 2 ** 64)
    if params['seed_g2'] == None:
        params['seed_g2'] = randint(0, 2 ** 64)
    if params['seed_g3'] == None:
        params['seed_g3'] = randint(0, 2 ** 64)

    if params['q']==None:
        ## for dihedral group, calculate for cyclic the power of two that gives
        ## error less than 2**-100, then q' = sqrt(2)*q for dihedral.
        n = params['n']
        ns = n/4 ##n^2
        smalln = int(sqrt(ns))
        d = int(ns/7)
        p =3
        q= nextprime(24*d*p)
        while ((q-1)%smalln!=0):
            q = nextprime(q)


        params['q'] = q

    if params['filename']==None:
        params['filename'] = str(params['n'])+'_'+str(params['q'])
    return params

def run_all(f, params):
    jobs = []

    file_tag = params['file_tag']
    params['group'] = 'bqtru'
    if params['dump']:
        create_file(file_tag, params['group'] , params['filename'], params['attack_type']) ##Create excel file to start saving the records
        #dump_seed(original_seed,params['group'],params['filename'])
    for t in range(params['trials']):
        params_  = copy.deepcopy(params)

        params['seed_f']  = random.randint(0, 2**64)
        params['seed_g0'] = random.randint(0, 2**64)
        params['seed_g1'] = random.randint(0, 2 ** 64)
        params['seed_g2'] = random.randint(0, 2 ** 64)
        params['seed_g3'] = random.randint(0, 2 ** 64)


        jobs.append(params_)
    if params['workers'] == 1:
        for job in jobs:
            res = f(copy.deepcopy(job))
    else:
        pool = Pool(params['workers'])
        pool.map(f, jobs)
