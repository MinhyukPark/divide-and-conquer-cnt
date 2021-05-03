import copy
import sys

def CN3(u,v):
    n = len(u)
    if len(v)!= n:
        print("Length doesn't match")
        return -1
    N = max(max(u),max(v))

    F = {}
    F[-1]={}
    for i in range(n):
        F[i] = {}


    # init
    for m in range(1,N+1):
        for p1 in range(N+1):
            for q1 in range(N+1):
                for p2 in range(N+1):
                    for q2 in range(N+1):
                        F[-1][m,p1,q1,p2,q2] = p1+q1+p2+q2

    previous_m = copy.deepcopy(F)
    for m in range(1,N+1):
        for p1 in range(N+1):
            for q1 in range(N+1):
                for p2 in range(N+1):
                    for q2 in range(N+1):
                        previous_m[-1][m,p1,q1,p2,q2] = None

    #step
    for i in range(n):
        #for m in range(1,N+1):
        for m in range(1,max(u[i],v[i])+1):
            max_p1 = N+1
            min_p1 = 0
            max_p2 = N+1
            min_p2 = 0
            if u[i]>0:
                max_p1 = m
                min_p1 = max(m-u[i],0)
            else: #u[i]=0
                min_p1 = m
            if v[i]>0:
                max_p2 = m
                min_p2 = max(m-v[i],0)
            else: #v[i]=0
                min_p2 = m
            for p1 in range(min_p1,max_p1):
                for p2 in range(min_p2,max_p2):
                    max_q1 = N+1
                    min_q1 = 0
                    max_q2 = N+1
                    min_q2 = 0
                    if u[i]>0:
                        min_q1 = u[i]-m+p1
                        max_q1 = min_q1 + 1
                    if v[i]>0:
                        min_q2 = v[i]-m+p2
                        max_q2 = min_q2 + 1
                    for q1 in range(min_q1,max_q1):
                        for q2 in range(min_q2,max_q2):
                            best = sys.maxsize
                            previous = None
                            for l,x1,y1,x2,y2 in F[i-1].keys():
                                current = F[i-1][l,x1,y1,x2,y2] + max(p1-x1,0) + max(q1-y1,0)+ max(p2-x2,0) + max(q2-y2,0)
                                if(current < best):
                                    best = current
                                    previous = (l,x1,y1,x2,y2)
                            F[i][m,p1,q1,p2,q2] = best
                            previous_m[i][m,p1,q1,p2,q2] = previous
    # print(F)

    d = sys.maxsize
    m_arr = []
    end_m = None
    for m,p1,q1,p2,q2 in F[n-1].keys():
        # print("p1: " + str(p1) + " p2: " + str(p2) + " m: " + str(m) + " q1: " + str(q1) + " q2: " + str(q2))
        if(F[n-1][m,p1,q1,p2,q2] < d):
            d = F[n-1][m,p1,q1,p2,q2]
            end_m = (m,p1,q1,p2,q2)

        # d = min(F[n-1][m,p1,q1,p2,q2],d)
    # ''' code for backtracking
    m_arr.append(end_m[0])
    # m,p1,q1,p2,q2 = end_m
    # next_m = previous_m[n-1][m,p1,q1,p2,q2]
    next_m = end_m
    # print("m_arr: " + str(m_arr))
    steps_arr = []
    for i in range(n-1, -1, -1):
        m,p1,q1,p2,q2 = next_m
        # print("next_m: " + str(next_m))
        steps_arr.append(next_m)
        current_entry = previous_m[i][m,p1,q1,p2,q2]
        # print("current_entry: " + str(current_entry))
        m_arr.append(current_entry[0])
        next_m = current_entry

    steps_arr = steps_arr[n-1::-1]
    # print(steps_arr)
    m_arr = m_arr[n-1::-1]
    cumulative_score = 0
    # print ("DP Distance = " + str(d))
    # print("final m = " + str(m_arr))
    # print("left = " + str(u))
    # print("right = " + str(v))
    for i in range(n):
        current_step = steps_arr[i]
        m,p1,q1,p2,q2 = current_step
        # print("p1: " + str(p1) + " p2: " + str(p2) + " m: " + str(m) + " q1: " + str(q1) + " q2: " + str(q2))
        # print("position " + str(i))
        # print("current task is " + str(m_arr[i]) + " branching to " + str(u[i]) + " and " + str(v[i]))
        # print("go from " + str(m_arr[i]) + " to " + str(u[i]) + " by " + str(p1) + " losses and " + str(q1) + " gains")
        # print("go from " + str(m_arr[i]) + " to " + str(v[i]) + " by " + str(p2) + " losses and " + str(q2) + " gains")
        current_score = F[i][m,p1,q1,p2,q2]
        # print("the cost of this move is " + str(current_score - cumulative_score))
        # print("the cost of all the moves so far is " + str(current_score))
        cumulative_score = current_score
    # '''

    return d,m_arr
