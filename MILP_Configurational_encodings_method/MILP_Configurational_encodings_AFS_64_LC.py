#coding=utf-8
import copy
import time
from gurobipy import *
from itertools import combinations,permutations


def mycallback(model,where):
    if where == GRB.Callback.MIP:
        objval = model.cbGet(GRB.Callback.MIP_OBJBST)
        if objval <= 16.5:
            model.terminate()


def gubi(file, round_num, rot, k, s, Seq):
    model = read(file)
    model.Params.PoolSearchMode = 2
    model.Params.PoolSolutions = 1

    if model.isMIP == 0:
        print('Model is not a MIP')
        exit(0)

    #model.optimize(mycallback)
    model.optimize()

    nSolution = model.SolCount
    print("Number of solutions found: " + str(nSolution))

    ####################################Output Result File Storage Path
    point1 = "./" + "64bit_" + str(s) + "_linear_" + str(round(round_num / 2)) + "Round_" + str(k) + ".txt"

    o = open(point1, 'a')
    o.write(str(s))
    o.write('\n')
    o.write(str(k))
    o.write(':')
    o.write('\n')
    o.write('\n')

    for e in range(nSolution):
        model.Params.SolutionNumber = e

        x_list = []
        for i in range(round_num + 1):
            x_list.append([0] * 64)

        y_list = []
        for i in range(round_num):
            y_list.append([0] * 32)


        p_list = []
        for i in range(round_num):
            p_list.append([0] * 33)


        print('--------------Solution:' + str(e) + '----------------' + '\n')

        # for v in model.getVars():
        #     if v.x == 1:
        #         print(v.varName, v.x)

        for v in model.getVars():
            name = v.varName
            value = int(v.x)

            name_split = name.split('_')
            # print(name)
            name_list = [name[0], name_split[0].replace(name[0], ''), name_split[1]]
            if value == 1 and name_list[0] == 'x':
                x_list[int(name_list[1])][int(name_list[2])] = value

            if value == 1 and name_list[0] == 'y':
                y_list[int(name_list[1])][int(name_list[2])] = value


            if value == 1 and name_list[0] == 'p':
                p_list[int(name_list[1])][int(name_list[2])] = value



        for i in range(round_num + 1):
            # print(f'x{i}:', str(x_list[i][::1]))
            o.write('x')
            o.write(str(i))
            o.write(':')
            o.write(str(x_list[i][::1]))
            o.write('\n')


        for i in range(round_num):
            # print(f'y{i}:', str(y_list[i][::1]))
            o.write('y')
            o.write(str(i))
            o.write(':')
            o.write(str(y_list[i][::1]))
            o.write('\n')



        for i in range(round_num):
            # print(f'p{i}:', str(p_list[i][::1]))
            o.write('p')
            o.write(str(i))
            o.write(':')
            o.write(str(p_list[i][::1]))
            o.write('\n')


    o.write('\n')


    if nSolution != 0:
        print('---------------------------------------')
        print('Obj:', model.ObjVal)
        print('---------------------------------------')
        m = model.ObjVal
        return m
    else:
        return -1



def rotl(A,a):          ######################circular left shift
    for i in range(a):
        A.insert(len(A),A[0])
        A.remove(A[0])
    return A


def rotr(A,a):             ######################circular right shift
    for i in range(a):
        A.insert(0,A.pop())
    return A


def branch_Constraint(a, b, c):     ########XOR constraint condition
    buf = ''
    buf = buf + ' + ' + str(a) + ' + ' + str(b) + ' - ' + str(c) + ' >= 0\n'
    buf = buf + ' + ' + str(a) + ' - ' + str(b) + ' + ' + str(c) + ' >= 0\n'
    buf = buf + ' - ' + str(a) + ' + ' + str(b) + ' + ' + str(c) + ' >= 0\n'
    buf = buf + ' - ' + str(a) + ' - ' + str(b) + ' - ' + str(c) + ' >= -2\n'
    return buf


def modulo_Constraint(A, B, C, P):  ###### modular addition constraint
    buf = ''
    for i in range(0, 32):
        buf = buf + str(P[i]) + ' - ' + str(C[i]) + ' - ' + str(A[i]) + ' + ' + str(B[i]) + ' + ' + str(P[i+1]) + ' >= 0\n'
        buf = buf + str(P[i]) + ' + ' + str(C[i]) + ' + ' + str(A[i]) + ' - ' + str(B[i]) + ' - ' + str(P[i+1]) + ' >= 0\n'
        buf = buf + str(P[i]) + ' + ' + str(C[i]) + ' - ' + str(A[i]) + ' - ' + str(B[i]) + ' + ' + str(P[i+1]) + ' >= 0\n'
        buf = buf + str(P[i]) + ' - ' + str(C[i]) + ' + ' + str(A[i]) + ' - ' + str(B[i]) + ' + ' + str(P[i+1]) + ' >= 0\n'

        buf = buf + str(P[i]) + ' + ' + str(C[i]) + ' - ' + str(A[i]) + ' + ' + str(B[i]) + ' - ' + str(P[i+1]) + ' >= 0\n'
        buf = buf + str(P[i]) + ' - ' + str(C[i]) + ' + ' + str(A[i]) + ' + ' + str(B[i]) + ' - ' + str(P[i+1]) + ' >= 0\n'
        buf = buf + str(C[i]) + ' - ' + str(P[i]) + ' + ' + str(A[i]) + ' + ' + str(B[i]) + ' + ' + str(P[i+1]) + ' >= 0\n'
        buf = buf + str(P[i]) + ' + ' + str(C[i]) + ' + ' + str(A[i]) + ' + ' + str(B[i]) + ' + ' + str(P[i+1]) + ' <= 4\n'
    return buf

def equal(A, B): #####linear mask directly propagates to the next round
    buf = ''
    for i in range(0,32):
        buf = buf + ' + ' + str(A[i]) + ' - ' + str(B[i]) + ' = 0\n'
    return buf


def linear(blocksize2, k1, round_num, rot, s1, Seq):

    #######Storage location for the output MILP model LP file.#######
    filename = "./" + "64bit_" + str(s1)+ "_linear_" + str(round(round_num/2))+ "Round_"+str(k1)+".lp"


    o = open(filename, 'w')

    X = [[] for i in range(round_num + 1)]
    P = [[] for i in range(round_num)]
    Y = [[] for i in range(round_num)]


    for i in range(round_num + 1):
        x0 = ['x' + str(i) + '_' + str(j) for j in range(blocksize2)]
        x1 = ['x' + str(i) + '_' + str(j) for j in range(blocksize2,   2*blocksize2)]
        X[i].append(x0)
        X[i].append(x1)

    for i in range(round_num):
        p0 = ['p' + str(i) + '_' + str(j) for j in range(blocksize2+1)]
        P[i].append(p0)

    for i in range(round_num):
        y0 = ['y' + str(i) + '_' + str(j) for j in range(blocksize2)]
        Y[i].append(y0)





    X_X = copy.deepcopy(X)


    o.write("Minimize\n")
    buf = ''
    for i in range(round_num):
        if s1[(i%8)] == 1:
            for j in range(1, blocksize2):
                buf = buf + "p" + str(i) + "_" + str(j) + " + "

    buf = buf[:-3]       
    o.write(buf)
    o.write('\n')
    o.write('\n')

    o.write("Subject to\n")
    for r in range(round_num):
        if s1[(r%8)] == 1:
            buf = ''  ################# Initial jump state is set to 0
            buf = buf + "p" + str(r) + '_0 = 0' + '\n'
            o.write(buf)

    for r in range(round_num):
        if s1[(r%8)] == 1:
            o.write(modulo_Constraint(X[r][0], Y[r][0], X[r+1][1], P[r][0]))

            rotl(Y[r][0], k1[(r % 8)])
            for i in range(blocksize2):
                o.write(branch_Constraint(Y[r][0][i], X[r+1][0][i], X[r][1][i]))

        if s1[(r%8)] == 0:
            rotl(X[r][0], k1[(r % 8)])
            for i in range(blocksize2):
                o.write(branch_Constraint(X[r][0][i], X[r+1][0][i], X[r][1][i]))

            rotr(X[r][0], k1[(r % 8)])
            o.write(equal(X[r+1][1], X[r][0]))




    buf = ''  ################# At least one bit of the initial input is active
    for j in range(0, 2*blocksize2):
        buf = buf + "x0_" + str(j)
        if j != (2*blocksize2-1):
            buf = buf + " + "
        if j == (2*blocksize2-1):
            buf = buf + " >= 1\n"
    o.write(buf)



    o.write("Bounds\n")
    o.write("Binary\n")

    buf = ''
    for i in range(round_num+1):
        for j in range(2*blocksize2):
            buf = buf + "x" + str(i) + "_" + str(j) + "\n"
    o.write(buf)

    buf = ''
    for i in range(round_num):
        if s1[(i%8)] == 1:
            for j in range(blocksize2 + 1):
                buf = buf + "p" + str(i) + "_" + str(j) + "\n"
    o.write(buf)

    buf = ''
    for i in range(round_num):
        if s1[(i%8)] == 1:
            for j in range(blocksize2):
                buf = buf + "y" + str(i) + "_" + str(j) + "\n"
    o.write(buf)



    o.write('End')
    o.close()

    result = gubi(filename, round_num, rot, k1, s1, Seq)
    return result


def set_value(a, value):
    buf = ''
    buf = buf + str(a) + " = " + str(value) + "\n"
    return buf



def main():
    blocksize2 = 32     ###S-box branch size
    N=4                 #### Number of Rounds of the S-box


    round_num = 2*N     ## Number of Rounds after the Decomposition of the S-box

################## s is the set of configuration encodings.
################## k is the set of rotation sequences.

    ##############Alzette(10101010)##################
    # s = [[1, 0, 1, 0, 1, 0, 1, 0]]
    # k = [[31, 24, 17, 17, 0, 31, 24, 16]]

    ##############AFS-64(11000011)##################
    s = [[1, 1, 0, 0, 0, 0, 1, 1]]
    k = [[1, 1, 16, 24, 31, 31, 8, 16],
         [1, 1, 16, 8, 31, 31, 24, 16],
         [17, 24, 1, 1, 16, 31, 24, 0],
         [17, 8, 1, 1, 16, 31, 8, 0]]






    for Seq in range(0, len(s)):
        for rot in range(0, len(k)):
            start = time.time()

            ####################################Output Result File Storage Path
            point = "./"+"64bit_"+ str(s[Seq])+ "_linear_"+ str(round(round_num/2)) + "Round_" + str(k[rot])+ ".txt"

            linear_flag = linear(blocksize2, k[rot], round_num, rot, s[Seq], Seq)
            point_file = open(point, "a")

            # if differ_flag >= -1:
            point_file.write(str(s[Seq]))
            point_file.write('\n')
            point_file.write(str(k[rot]))
            linear_flag1 = repr(linear_flag)

            point_file.write("linear correlation : " + linear_flag1)
            point_file.write('\n')
            point_file.write('\n')
            point_file.close()


            end = time.time()
            time1 = end - start
            point_file = open(point, "a")
            point_file.write("total_time : " + str(time1) + 's')
            point_file.write('\n')
            point_file.write('\n')
            point_file.close()


if __name__ == '__main__':
    main()
