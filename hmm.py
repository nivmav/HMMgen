import os
import sys
import random
import math
import pickle


DELTA = 0.001
class HMM:
    N = 0
    M = 0
    A = []
    B = []
    pi =[]
    alpha = []
    beta = []
    gamma = []
    xi = []
    di = []
    seqMap = []
    initpprob = 0
    pprob = 0
    C = []
    T = 0
    L = 0
    Obs = []

    def Input(self, file1, file2):
        [pid1, seq1] = zip(*[line.split(' ') for line in open(file1)])
        [pid2, seq2] = zip(*[line.split(' ') for line in open(file2)])

        pid = list(pid1) + list(pid2)
        seq = list(seq1) + list(seq2)
        seq = map(lambda s: s.strip(), seq)
        self.T = len(seq)
        seen = set()
        seqUnique = [ x for x in seq if str( x ) not in seen and not seen.add( str( x ) )]


        self.di = dict(zip(seqUnique, xrange(len(seqUnique))))
        self.seqMap = map(lambda s: self.di[s], seq)

        return pid, seq, seqUnique

    def getProcSeq(self):

        seq = self.Obs

        seen = set()
        seqUnique = [ x for x in seq if str( x ) not in seen and not seen.add( str( x ) )]

        self.N = 22
        self.M = 22

        self.di = dict(zip(seqUnique, xrange(len(seqUnique))))
        self.seqMap = map(lambda s: self.di[s], seq)

        self.T = len(self.seqMap)
        self.C = [0.0 for i in range(self.T)]
        self.alpha = [ [ 0.0 for j in range(self.N) ] for i in range(self.T) ]
        self.HMMforward()

    def Data(self, name):
        from itertools import groupby
        things = [map(lambda s:s.strip(), line.split(' ')) for line in open(name)]

        da = {}
        for key, group in groupby(things, lambda x: x[0]):
            lis = []
            for thing in group:
                lis.append(thing[1])
            da[key] = lis
        return da

    def CoreBase(self, train1, train2, test):
        import operator
        files = [train1, train2, test]
        base = {}
        count = 0
        for f in files:
            print files[count]
            count += 1
            data = self.Data(f)
            for key in data.keys():
                self.Obs = data[key]
                self.getProcSeq()
                base[key] = self.pprob
                print "likelihood %s is %d" % (key, self.pprob)
        self.sort_result = sorted(base.iteritems(), key = operator.itemgetter(1))


    def randMat(self, length):
        m = [ [ 0 for i in range(length) ] for j in range(length) ]
        for i in range(length):
            count = 0
            for j in range(length):
                m[i][j] = random.randrange(1,10)
                m[i][j] = float(m[i][j])
                count = count + m[i][j]
            for j in range(length):
                m[i][j] = m[i][j]/count

        return m
    def randArray(self, length):
        m = []
        count = 0
        for i in range(0, length):
            a = random.randrange(1, 10)
            m.append(float(a))
            count = count + a 
        for i in range(0, length):
            m[i] = m[i]/count


        return m

    def HMMforward(self):
        
        #Initialization
        self.C[0] = 0.0
        for i in range(0, self.N):
            self.alpha[0][i] = self.pi[i] * (self.B[i][self.seqMap[0]]);
            self.C[0] += self.alpha[0][i]

        for i in range(0, self.N):
            self.alpha[0][i] /= self.C[0]

        #Induction
        for t in range(0, self.T - 1):
            self.C[t+1] = 0.0
            for j in range(0,self.N):
                count = 0.0
                for i in range(0,self.N):
                    count = count + self.alpha[t][i] * self.A[i][j]
                self.alpha[t+1][j] = count * self.B[j][self.seqMap[t+1]]
                self.C[t+1] += self.alpha[t+1][j];

            for j in range(0, self.N):
                self.alpha[t+1][j] /= self.C[t+1]
        

        #Termination

        self.pprob = 0.0

        for t in range(0, self.T):
            self.pprob += math.log(self.C[t])
        

    def HMMbackward(self):
    
        #Initialization

        for i in range(0, self.N):
            self.beta[self.T-1][i] = 1.0/self.C[self.T-1]

        #Induction
        for t in xrange(self.T-2, -1, -1):
            for i in range(0, self.N):
                count = 0.0
                for j in range(0, self.N):
                    count += self.A[i][j] * (self.B[j][self.seqMap[t+1]])*self.beta[t+1][j]

                self.beta[t][i] = count/self.C[t];
    
    def ComputeG(self):
        
        for t in range(0, self.T):
            denominator = 0.0
            for j in range(0, self.N):
                self.gamma[t][j] = self.alpha[t][j]*(self.beta[t][j])
                denominator += self.gamma[t][j]



            for i in range(0, self.N):
                self.gamma[t][i] = self.gamma[t][i]/denominator


    def ComputeXi(self):
        for t in range(0, self.T-1):
            count = 0.0
            for i in range(0, self.N):
                for j in range(0, self.N):
                    self.xi[t][i][j] = self.alpha[t][i] *(self.beta[t+1][j]) * (self.A[i][j])*(self.B[j][self.seqMap[t+1]])
                    count += self.xi[t][i][j]

            for i in range(0, self.N):
                for j in range(0, self.N):
                    self.xi[t][i][j] /= count





    def BaumWelch(self):


        self.C = [0.0 for i in range(self.T)]
        self.alpha = [ [ 0.0 for j in range(self.N) ] for i in range(self.T) ]

        self.beta = [ [ 0.0 for j in range(self.N) ] for i in range(self.T) ]


        deltaprev = 10e-70
        self.gamma = [ [0.0  for j in range(self.N)] for i in range(self.T)]
        self.xi = []
        for t in range(self.T):
            dxi = [[0.0 for j in range(0, self.N)] for i in range(0, self.N)]
            self.xi.append(dxi)
        scale = [0.0 for i in range(0, self.T)]
        self.HMMforward()
        loginit = self.pprob
        self.initpprob = self.pprob
        self.HMMbackward()
        self.ComputeG()
        self.ComputeXi()
        self.L = 0
        while True:
            for i in range(0, self.N):
                self.pi[i] == .001 + .999*self.gamma[0][i];

            for i in range(0, self.N):
                denominatorA = 0.0
                for t in range(0, self.T-1):
                    denominatorA += self.gamma[t][i]

                for j in range(0, self.N):
                    numeratorA = 0.0
                    for t in range(0, self.T-1):
                        numeratorA += self.xi[t][i][j]
                    self.A[i][j] = .001 + .999*(numeratorA/denominatorA)


                denominatorB = denominatorA + self.gamma[self.T-1][i]
                for k in range(0, self.M):
                    numeratorB = 0.0
                    for t in range(0, self.T):
                        if (self.seqMap[t] == k):
                            numeratorB += self.gamma[t][i]

                    self.B[i][k] = .001 + .999*numeratorB/denominatorB

            
            self.HMMforward()
            self.HMMbackward()
            self.ComputeG()
            self.ComputeXi()
                
            print "In a loop pprob  %d" % self.pprob
            print "In a  loop loginit %d" % loginit
            delta = self.pprob - loginit
            loginit = self.pprob
            print delta
            
            self.L += 1
            
            if(delta < DELTA):
                break
        print self.L
        print self.pprob
        print self.initpprob
    



    def Process(self):

        self.pid, self.seq, self.seqUnique = self.Input("list1.txt", "list2.txt");
        self.N = len(self.seqUnique)
        self.M = self.N

        self.A = self.randMat(self.N)
        self.B = self.randMat(self.N)
        self.pi = self.randArray(self.N)
        self.BaumWelch()

    def Unpickle(self):
        self.A = pickle.load(open("a.lambda", "rb"))
        self.B = pickle.load(open("b.lambda", "rb"))
        self.pi = pickle.load(open("pi.lambda", "rb"))



def main():
    
    hmmnew = HMM()
    hmmnew.Process()

    #hmmnew.Unpickle()
    hmmnew.CoreBase("list1.txt", "list2.txt", "test.txt")
    #print pid
    #print seq





if __name__ == "__main__":
    sys.exit(main())
