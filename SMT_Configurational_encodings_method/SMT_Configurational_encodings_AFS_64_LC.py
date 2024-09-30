import subprocess
import re
import sys
import time


def AFS_P(round, cvc_path1, prob, s, k, s_l, k_l):
    ROUNDS = round
    fp = open(cvc_path1, "w")

    for r in range(ROUNDS+1):

      
        fp.write("x_{}: BITVECTOR(32) ;".format(r))  # Left input linear characteristic for the r-th round: x_r;
        fp.write("\n")
       
        fp.write("y_{}: BITVECTOR(32) ;".format(r))  # Right input linear characteristic for the r-th round: y_r;
        fp.write("\n")


   
    for r in range(ROUNDS):
        if s[r%s_l] == 1:
            fp.write("p_{}: BITVECTOR(32) ;".format(r))  # Correlation auxiliary variable for modular addition;
            fp.write("\n")

            fp.write("z_{}: BITVECTOR(32) ;".format(r))  # Linear characteristic before the rotation operation in the r-th round: z_r
            fp.write("\n")
  
    fp.write("total_probability: BITVECTOR(32) ;")  # Total target Correlation;
    fp.write("\n")
    fp.write("ASSERT( NOT( (x_0 @ y_0) = 0hex0000000000000000) );\n") # Input Linear characteristic are not all zero;


    for r in range(ROUNDS):


        if s[r%s_l] == 1:
            if k[r%k_l] != 0:

                fp.write("ASSERT( x_{}=BVXOR(y_{},z_{}) ) ;\n".format(r+1,r,r))  ## Linear propagation of branch operation： x_r+1 = ( y_r xor z_r )

                ################################ p_r = (M_n)^T(y_r+1,x_r,z_r)
                fp.write("ASSERT( p_" + str(r) + " = (0bin0)@(BVXOR(BVXOR( x_{},(z_{}[{}:0])@(z_{}[31:{}]) ), y_{})[{}:{}])".format(r, r, k[r % k_l]-1, r, k[r % k_l], r+1, 31, 31))
                for i in range(2,32): #1~n-1,   n=32
                    fp.write("@(")
                    for j in range(2, i+1):
                        fp.write("BVXOR(")
                    fp.write("BVXOR(BVXOR( x_{},(z_{}[{}:0])@(z_{}[31:{}]) ), y_{})[{}:{}]".format(r,r,k[r % k_l]-1,r,k[r % k_l],r+1,31,31))
                    for j in range(2, i + 1):
                        fp.write(",BVXOR(BVXOR( x_{},(z_{}[{}:0])@(z_{}[31:{}]) ), y_{})[{}:{}] )".format(r,r,k[r % k_l]-1,r,k[r % k_l],r+1,32-j,32-j))
                    fp.write(")")
                fp.write(") ;\n")
                ##############################

                fp.write("ASSERT( BVXOR(y_{},x_{}) = ( BVXOR(y_{},x_{})&p_{}) ) ;\n".format(r+1,r,r+1,r,r)) ### Condition for the effective linear feature propagation in modular arithmetic operations： ( y_r+1 xor x_r ) is less than in the partial order p_r ====> ( y_r+1 xor x_r ) =  ( y_r+1 xor x_r ) & p_r
                fp.write("ASSERT( BVXOR( y_{},(z_{}[{}:0])@(z_{}[31:{}]) ) = ( BVXOR( y_{},(z_{}[{}:0])@(z_{}[31:{}]) )&p_{}) ) ;\n".format(r + 1, r,k[r % k_l]-1,r,k[r % k_l],r+1,r,k[r % k_l]-1,r,k[r % k_l],r) )


            if k[r%k_l] == 0:

                fp.write("ASSERT( x_{}=BVXOR(y_{},z_{}) ) ;\n".format(r+1,r,r))  ##linear propagation of branch operations： x_r+1 = ( y_r xor z_r )

                ################################ p_r = (M_n)^T(y_r+1,x_r,z_r)
                fp.write("ASSERT( p_" + str(r) + " = (0bin0)@(BVXOR(BVXOR( x_{}, z_{} ), y_{})[{}:{}])".format(r, r, r+1, 31, 31))
                for i in range(2,32): #1~n-1,   n=32
                    fp.write("@(")
                    for j in range(2, i+1):
                        fp.write("BVXOR(")
                    fp.write("BVXOR(BVXOR( x_{}, z_{} ), y_{})[{}:{}]".format(r,r,r+1,31,31))
                    for j in range(2, i + 1):
                        fp.write(",BVXOR(BVXOR( x_{}, z_{} ), y_{})[{}:{}] )".format(r,r,r+1,32-j,32-j))
                    fp.write(")")
                fp.write(") ;\n")
                ##############################

                fp.write("ASSERT( BVXOR(y_{},x_{}) = ( BVXOR(y_{},x_{})&p_{}) ) ;\n".format(r+1,r,r+1,r,r)) ### Condition for the effective linear feature propagation in modular arithmetic operations： ( y_r+1 xor x_r ) is less than in the partial order p_r ====> ( y_r+1 xor x_r ) =  ( y_r+1 xor x_r ) & p_r
                fp.write("ASSERT( BVXOR( y_{}, z_{} ) = ( BVXOR( y_{}, z_{} )&p_{}) ) ;\n".format(r + 1, r,r+1,r,r) )


        if s[r%s_l] == 0:

            if k[r%k_l] != 0:
                fp.write("ASSERT(x_{} = y_{}) ;\n".format(r, r+1)) ##linear propagation of XOR operations： x_r = y_(r+1)

                fp.write("ASSERT( x_{} = BVXOR(y_{}, (x_{}[{}:0])@(x_{}[31:{}]) ) ) ;\n".format(r+1, r, r, 31-k[r % k_l], r, 32-k[r % k_l]))  ##linear propagation of branch operations: x_r+1 = y_r xor (x_r<<k[r])

            if k[r%k_l] == 0:
                fp.write("ASSERT(x_{} = y_{}) ;\n".format(r, r+1)) ##linear propagation of XOR operations： x_r = y_(r+1)

                fp.write("ASSERT( x_{} = BVXOR(y_{}, x_{} ) ) ;\n".format(r+1, r, r))  ##linear propagation of branch operations: x_r+1 = y_r xor (x_r<<k[r])


        
    



    fp.write("ASSERT( total_probability = BVPLUS(32")
    for r in range(ROUNDS):
        if s[r%s_l] == 1:
            for i in range(0, 31): ####Addition of correlation-assisting bits from the (n-2)th position to the least significant bit at position 0.
                fp.write(", (0bin0000000000000000000000000000000@(p_{}[{}:{}]))".format(r, i, i))
    fp.write(" ) ) ; \n")  # 概率相加

    fp.write("ASSERT( total_probability = 0hex" + f"{prob:08x}" + " ) ;\n")  # Target correlation constraint where 'prob' is the exponent of the correlation.
    fp.write("QUERY(FALSE) ;\n")
    fp.write("COUNTEREXAMPLE ;\n")
    fp.close()

def run_stp(r, cvc_path, result_path, s, k, s_l, k_l):
    global result, prob,low,upper
    start_time = time.time()
    round = r

    cvc_path1 = cvc_path

    
    try:

        stp_parameters = ["stp", cvc_path]


        for prob in range(low, upper+1, 1):
        
            AFS_P(round, cvc_path1, prob, s, k, s_l, k_l)
            result = subprocess.check_output(stp_parameters, universal_newlines=True, stderr=subprocess.STDOUT)
            if "Valid." in result:
                print("Not Find: 2^-",prob)

            else:
                print("Find: 2^-", prob)
                print()
                break

        if prob>=0:
            with open(result_path, 'w') as f:
                f.write("Linear Correlation:2^-" + str(prob) + '\n')
                f.write("Linear Characters：\n")

                pattern1 = r'x_\d+ = 0x[\da-fA-F]+'
                matches1 = re.findall(pattern1, result)
                sorted_matches1 = sorted(matches1, key=lambda x: int(x.split('_')[1].split(' ')[0]))
                for match1 in sorted_matches1:
                    f.write(match1 + '\n')

                pattern2 = r'y_\d+ = 0x[\da-fA-F]+'    
                matches2 = re.findall(pattern2, result)
                sorted_matches2 = sorted(matches2, key=lambda y: int(y.split('_')[1].split(' ')[0]))
                for match2 in sorted_matches2:
                    f.write(match2 + '\n')
                
                pattern3 = r'p_\d+ = 0x[\da-fA-F]+'    
                matches3 = re.findall(pattern3, result)
                sorted_matches3 = sorted(matches3, key=lambda p: int(p.split('_')[1].split(' ')[0]))
                for match3 in sorted_matches3:
                    f.write(match3 + '\n')


                f.write("\n")
                end_time = time.time()
                f.write("Running Time：" + str(end_time - start_time) + "秒" + '\n')

                f.close()


    except subprocess.CalledProcessError as e:
        print(f"STP error: {e}")
        return False



low=0   # Target lower bound
#upper = low  # Fixed single target
upper =9999 # Target upper bound
def main():

    ### s_pool is the set of configuration encodings.
    ### k_pool is the set of rotation sequences.

    ###############AFS-64(11000011)###################
    s_pool = [[1, 1, 0, 0, 0, 0, 1, 1]]
    k_pool = [[1, 1, 16, 24, 31, 31, 8, 16], [1, 1, 16, 8, 31, 31, 24, 16], [17, 24, 1, 1, 16, 31, 24, 0], [17, 8, 1, 1, 16, 31, 8, 0]]

    # ###############Alzette(10101010)###################
    # s_pool = [[1, 0, 1, 0, 1, 0, 1, 0]]
    # k_pool = [[31,24,17,17,0,31,24,16]]

    rd=4       ### S-box Number of rounds

    for s_i in range(len(s_pool)):
        for k_i in range(len(k_pool)):

            s = s_pool[s_i]
            s_l = len(s)
            # s_str = ','.join(map(str, s))
            k = k_pool[k_i]
            k_l = len(k)
            # k_str = ','.join(map(str,k))

            PATHTEMP1 = "./" + "64bit_LC_s_" + str(s) + "k_" + str(k) + "r_" + str(rd) + "low_" + str(low) + ".cvc"     ###Output/Input SMT Model CVC File Storage Path
            PATHTEMP2 = "./" + "64bit_LC_s=" + str(s) + ";k=" + str(k) + ";r=" + str(rd) + ";low=" + str(low) + ".txt"  ###Output Result File Storage Path

            run_stp(2*rd, PATHTEMP1, PATHTEMP2, s, k, s_l, k_l) # (Number of rounds after structural decomposition, path to store model .cvc files, path to store output result files)

if __name__ == '__main__':
    main()
