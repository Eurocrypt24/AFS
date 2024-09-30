import subprocess
import re
import sys
import time
import ast

def AFS_P(round, cvc_path1, prob, s, k, s_l, k_l):
    ROUNDS = round
    fp = open(cvc_path1, "w")

    for r in range(ROUNDS+1):
     
      
        fp.write("x_{}: BITVECTOR(16) ;".format(r))  # Left input difference for the r-th round: x_r;
        fp.write("\n")
       
        fp.write("y_{}: BITVECTOR(16) ;".format(r))  # Right input difference for the r-th round: y_r;
        fp.write("\n")
   
    for r in range(ROUNDS):
        if s[r%s_l] == 1:
            fp.write("p_{}: BITVECTOR(16) ;".format(r))  # Probability auxiliary variable for modular addition;
            fp.write("\n")
  
    fp.write("total_probability: BITVECTOR(16) ;")  # Total target probability;
    fp.write("\n")
    fp.write("ASSERT( NOT( (x_0 @ y_0) = 0hex00000000) );")  # Input differences are not all zero;


    for r in range(ROUNDS):
    

        if s[r%s_l] == 1:
            if k[r%k_l] != 0:
               
                fp.write("ASSERT( p_{} = BVXOR( (~(x_{} << 1)[15:0]), (((y_{}[{}:0])@(y_{}[15:{}])) << 1)[15:0] ) & BVXOR( (~(x_{} << 1)[15:0]), (y_{} << 1)[15:0] ) );\n".format(r, r, r, k[r%k_l]-1, r, k[r%k_l], r, r+1))  #eq(a<<1,b<<1,c<<1) = (~(a<<1)xor(b<<1))&(~(a<<1)xor(c<<1))
                fp.write("ASSERT( p_{} & BVXOR(BVXOR(BVXOR(x_{}, (y_{}[{}:0])@(y_{}[15:{}]) ), y_{}), (((y_{}[{}:0])@(y_{}[15:{}])) << 1)[15:0] ) = 0hex0000 ) ;\n".format(r, r, r, k[r%k_l]-1, r, k[r%k_l], r+1, r, k[r%k_l]-1, r, k[r%k_l])) # Conditions for effective differential propagation in modular addition operations
            
            if k[r%k_l] == 0:
                
                fp.write("ASSERT( p_{} = BVXOR( ~(x_{} << 1)[15:0], (y_{} << 1)[15:0] ) & BVXOR( ~(x_{} << 1)[15:0], (y_{} << 1)[15:0] ) );\n".format(r, r, r, r, r+1))  #eq(a<<1,b<<1,c<<1) = (~(a<<1)xor(b<<1))&(~(a<<1)xor(c<<1))
                fp.write("ASSERT( p_{} & BVXOR(BVXOR(BVXOR(x_{}, y_{}), y_{}), (y_{} << 1)[15:0] ) = 0hex0000 ) ;\n".format(r, r, r, r+1, r))   # Conditions for effective differential propagation in modular addition operations

            fp.write("ASSERT( x_{} = y_{} ) ;\n".format(r+1, r)) # Differential propagation for branch operation: a=b=c

        if s[r%s_l] == 0:
            if k[r%k_l] != 0:
                fp.write("ASSERT( y_{} = BVXOR(x_{}, (y_{}[{}:0])@(y_{}[15:{}]) ) ) ;\n".format(r+1, r, r, k[r%k_l]-1, r, k[r%k_l])) # Differential propagation in XOR operation: a xor b = c
                fp.write("ASSERT( x_{} = y_{} ) ;\n".format(r+1, r))  # Differential propagation for branch operation: a=b=c

            if k[r%k_l] == 0:
                fp.write("ASSERT( y_{} = BVXOR(x_{}, y_{}) ) ;\n".format(r+1, r, r)) # Differential propagation in XOR operation: a xor b = c
                fp.write("ASSERT( x_{} = y_{} ) ;\n".format(r+1, r))  # Differential propagation for branch operation: a=b=c


        
    



    fp.write("ASSERT( total_probability = BVPLUS(16")
    for r in range(ROUNDS):
        if s[r%s_l] == 1:
            for i in range(1, 16): #p_{}=eq(a<<1, b<<1, c<<1),Therefore, the most significant bit has shifted from position n-1 to the least significant position 0. By calculating the number of bits that are unequal from position n-1 down to position 1, we can determine the exponent for the probability calculation.
                fp.write(", (0bin000000000000000@(~p_{}[{}:{}]))".format(r, i, i))
    fp.write(" ) ) ; \n")

    fp.write("ASSERT( total_probability = 0hex" + f"{prob:04x}" + " ) ;\n")  # Target probability constraint where 'prob' is the exponent of the probability
    fp.write("QUERY(FALSE) ;")
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
                f.write("Differential Probability:2^-" + str(prob) + '\n')
                f.write("Differential Characters：\n")

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
                f.write("Running Time：" + str(end_time - start_time) + "s" + '\n')

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



    ###############AFS-32(11000011)###################
    s_pool = [[1, 1, 0, 0, 0, 0, 1, 1]]
    k_pool = [[9,8,2,1,7,15,14,0],[0,1,9,15,14,8,2,7]]

    # ###############AFS-32(10100011)###################
    # s_pool = [[1, 0, 1, 0, 0, 0, 1, 1]]
    # k_pool = [[0,7,14,15,8,1,2,9],[8,0,2,1,9,14,7,15]]

    rd=4  ### S-box Number of rounds

    for s_i in range(len(s_pool)):
        for k_i in range(len(k_pool)):
            s=s_pool[s_i]
            s_l=len(s)

            k=k_pool[k_i]
            k_l=len(k)


            PATHTEMP1 = "./" + "32bit_DP_s_" + str(s) + "k_" + str(k) + "r_" + str(rd) + "low_" + str(low) + ".cvc"  ###Output/Input SMT Model CVC File Storage Path
            PATHTEMP2 = "./" + "32bit_DP_s=" + str(s) + ";k=" + str(k) + ";r=" + str(rd) + ";low=" + str(low) + ".txt"  ###Output Result File Storage Path

            run_stp(2*rd, PATHTEMP1, PATHTEMP2, s, k, s_l, k_l) # (Number of rounds after structural decomposition, path to store model .cvc files, path to store output result files)

if __name__ == '__main__':
    main()
