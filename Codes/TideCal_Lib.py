# This is the library, including the packages, functions and PowerSystem Class, which are used in TideCal_Final.py.
from math import atan as atan
from math import pi as pi
from math import sin as sin
from math import cos as cos
from scipy import linalg
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
matplotlib.rcParams['font.sans-serif']=['Times New Roman']
epsilon=1e-15
class PowerSystem:
    def layout(self):#节点拓扑
        plt.figure(figsize=(8, 8))
        n=self.n
        G = nx.Graph()
        for i in range(1,n+1):
            G.add_node(i)
        for edge in self.edge:
            G.add_edge(edge[0],edge[1])
        pos=nx.spectral_layout(G)
        nx.draw(G,pos)
        node_label={}
        for node in G.nodes:
            node_label[node]=node
        nx.draw_networkx_labels(G,pos,labels=node_label)
        plt.savefig('1 layout.png',dpi=1000)
        plt.show()

    def layout2(self):#潮流分布
        plt.figure(figsize=(10,10))
        n = self.n
        G = nx.DiGraph()
        for i in range(1, n + 1):
            G.add_node(i)
        node_label = {}
        #edge_label = {}
        for node in G.nodes:
            node_label[node] = str(node)+'\nU='+printcomplex2(self.U[node-1])+'\nS='+printcomplex2(self.S[node-1])
        for edge in self.edge:
            key1 = str(edge[0]) + ',' + str(edge[1])
            key2 = str(edge[1]) + ',' + str(edge[0])
            if key1 in self.Sbranch:
                key=key1
            else:
                key=key2
            u,v=map(int,key.split(','))
            G.add_edge(u,v,S=self.Sbranch[key])
        edge_label=nx.get_edge_attributes(G,'S')
        pos = nx.spectral_layout(G)
        nx.draw(G, pos)
        nx.draw_networkx_labels(G, pos, labels=node_label)
        nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_label)
        plt.savefig('2 layout.png', dpi=1000)
        plt.show()

    def init_fromsile(self,path1,path2,path3):#读取文件数据，变量初始化，生成导纳矩阵
        #网络构建
        #文件读取
        # print("输入节点、边的数量和变压器数量")

        branch=[]
        with open(path1) as f:
            line = f.readline()
            while line:
                branch.append(list(map(eval,line.split())))
                line=f.readline()
            n=0
            for branchi in branch:
                if max(branchi[0],branchi[1])>n:
                    n=max(branchi[0],branchi[1])
            self.n=n
            self.edge = []
            self.y0 = [0] * n
            self.y = [[0] * n for _ in range(n)]
            self.U = [0] * n
            self.S = [0] * n
            # print("输入边和边的导纳(如1 2 5-15j)，标号请从1开始")
            for branchi in branch:
                start,end,R,L,S=branchi[0],branchi[1],branchi[2],branchi[3],branchi[4]
                self.edge.append([start,end])
                self.y[start-1][end-1]=1/(R+1j*L)
                self.y[end-1][start-1]=1/(R+1j*L)
                self.y0[start - 1] += 1j*S / 2
                self.y0[end - 1] += 1j*S / 2

        self.inputNode=[]#发电机接入的节点
        with open(path2) as f:
            line = f.readline()
            while line:
                line=list(map(eval,line.split()))
                self.inputNode.append(line[1])
                line = f.readline()

        with open(path3) as f:
            line = f.readline()
            while line:
                line = list(line.split())
                index=int(line[0])-1
                if line[1]!='?':
                    A,theta=map(eval,[line[1],line[2]])
                    self.U[index]=A*(cos(theta/180*pi)+1j*sin(theta/180*pi))
                else:
                    if int(line[0]) in self.inputNode:
                        P, Q = map(eval, [line[3], line[4]])
                        self.S[index] = P + 1j * (Q+epsilon)
                    else:
                        P, Q = map(eval, [line[5], line[6]])
                        self.S[index] =- (P + 1j * (Q+epsilon))
                line = f.readline()
            self.YB = [[0] * self.n for _ in range(self.n)]
            self.YB_Cal(self.y)
            self.dP = [0] * n
            self.dQ = [0] * n
            self.p = [0] * n
            self.q = [0] * n
            self.a = [[0] * n for _ in range(n)]
            self.b = [[0] * n for _ in range(n)]
            self.Jaccob = [[0] * (2 * n - 2) for _ in range(2 * n - 2)]
            self.H = [[0] * (n - 1) for _ in range(n - 1)]
            self.N = [[0] * (n - 1) for _ in range(n - 1)]
            self.J = [[0] * (n - 1) for _ in range(n - 1)]
            self.L = [[0] * (n - 1) for _ in range(n - 1)]
            self.cnt = 0

    def YB_Cal(self,y):
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    self.YB[i][j] = sum(y[i])+self.y0[i]
                else:
                    self.YB[i][j] = -y[i][j]

    def Jacood_Cal(self):
        # 雅可比元素计算
        n=self.n
        for i in range(2, n+1):
            Ii = (self.P(i) - 1j * self.Q(i)) / (self.U[i - 1].conjugate())
            self.a[i - 1][i - 1] = Re(Ii)
            self.b[i - 1][i - 1] = Im(Ii)

        for i in range(2, n+1):
            for j in range(2, n+1):
                self.H[i - 2][j - 2] = -self.B(i, j) * self.e(i) + self.G(i, j) * self.f(i) + self.b[i - 1][j - 1]
                self.N[i - 2][j - 2] = +self.G(i, j) * self.e(i) + self.B(i, j) * self.f(i) + self.a[i - 1][j - 1]
                self.J[i - 2][j - 2] = -self.G(i, j) * self.e(i) - self.B(i, j) * self.f(i) + self.a[i - 1][j - 1]
                self.L[i - 2][j - 2] = -self.B(i, j) * self.e(i) + self.G(i, j) * self.f(i) - self.b[i - 1][j - 1]
                # print("H {} {}".format(i, j))
                # print("%.4f" % H[i - 2][j - 2])
                # print("N {} {}".format(i, j))
                # print("%.4f" % N[i - 2][j - 2])
                # print("J {} {}".format(i, j))
                # print("%.4f" % J[i - 2][j - 2])
                # print("L {} {}".format(i, j))
                # print("%.4f" % L[i - 2][j - 2])
        #print(np.array(self.H))
        for i in range(2*n-2):
            for j in range(2*n-2):
                if i % 2 == 0 and j % 2 == 0:
                    self.Jaccob[i][j] = self.H[i // 2][j // 2]
                elif i % 2 == 0 and j % 2 != 0:
                    self.Jaccob[i][j] = self.N[i // 2][j // 2]
                elif i % 2 != 0 and j % 2 == 0:
                    self.Jaccob[i][j] = self.J[i // 2][j // 2]
                else:
                    self.Jaccob[i][j] = self.L[i // 2][j // 2]

    def renew(self):#further...
        ...

    def initialize(self):#初始化
        n=self.n
        #print(np.array(self.y))
        #print(np.array([[0,5-15j,1.25-3.75j,0,0],[5-15j,0,1.667-5j,1.667-5j,2.5-7.5j],[1.25-3.75j,1.667-5j,0,10-30j,0],[0,1.667-5j,10-30j,0,1.25-3.75j],[0,2.5-7.5j,0,1.25-3.75j,0]]))
        #print(np.array(self.YB))
        cnt=self.cnt
        for i in range(self.n):
            if self.U[i]!=0:
                index=i
                break

        for i in range(self.n):
            if i!=index:
                self.U[i]=self.U[0]
        print("初始化...")
        # 不平衡量计算
        for i in range(1, n+1):
            self.p[i - 1] = self.P(i)
            #print("P{}({}):".format(i, cnt))
            #print("%.4f" % p[i - 1])
            self.q[i - 1] = self.Q(i)
            #print("Q{}({}):".format(i, cnt))
            #print("%.4f" % q[i - 1])

        for i in range(2, n+1):
            self.dP[i - 1] = Re(self.S[i - 1]) - self.P(i)
            self.dQ[i - 1] = Im(self.S[i - 1]) - self.Q(i)
            # print("dP{}({}):".format(i, cnt))
            # print("%.4f" % self.dP[i - 1])
            # print("dQ{}({}):".format(i, cnt))
            # print("%.4f" % self.dQ[i - 1])
        self.Jacood_Cal()
        #print(np.array(self.Jaccob))
        # 求节点电压
        self.Y = []
        for i in range(1, n):
            self.Y.append(self.dP[i])
            self.Y.append(self.dQ[i])

        self.de_f = np.dot(linalg.inv(self.Jaccob), np.array(self.Y))
        for i in range(2*n-2):
            if i % 2 == 0:
                self.U[i // 2 + 1] += 1j * self.de_f[i]
            else:
                self.U[i // 2 + 1] += self.de_f[i]
        print("雅可比矩阵")
        print(self.Jaccob)
        print("节点功率不平衡量")
        print(self.Y)
        print("节点电压修正量")
        print(self.de_f)
        print("节点电压新值")
        print(self.U)
        print("")

    def iterate(self):#迭代，直至误差小于阈值
        n=self.n
        while norm(self.de_f) > 1e-5:
            self.cnt += 1
            print("Times:%d" % self.cnt)

            # 不平衡量计算
            for i in range(1, n+1):
                self.p[i - 1] = self.P(i)
                # print("P{}({}):".format(i, self.cnt))
                # print("%.10f" % self.p[i - 1])
                self.q[i - 1] = self.Q(i)
                # print("Q{}({}):".format(i, cnt))
                # print("%.4f" % q[i - 1])

            for i in range(2, n+1):
                self.dP[i - 1] = Re(self.S[i - 1]) - self.P(i)
                self.dQ[i - 1] = Im(self.S[i - 1]) - self.Q(i)
                # print("dP{}({}):".format(i, self.cnt))
                # print("%.4f" % self.dP[i - 1])
                # print("dQ{}({}):".format(i, self.cnt))
                # print("%.4f" % self.dQ[i - 1])

            self.Jacood_Cal()
            # 求节点电压
            # print(linalg.inv(Jaccob))
            # Y = []
            for i in range(1, n):
                self.Y[i * 2 - 2] = self.dP[i]
                self.Y[i * 2 - 1] = self.dQ[i]
            #print("Y", self.Y)
            self.de_f = np.dot(linalg.inv(self.Jaccob), np.array(self.Y))
            #J_inv=linalg.inv(self.Jaccob)
            #self.de_f=[sum([J_inv[i][j]*self.Y[j] for j in range(2*n-2)]) for i in range(2*n-2)]
            for i in range(2*n-2):
                if i % 2 == 0:
                    self.U[i // 2 + 1] += 1j * self.de_f[i]
                else:
                    self.U[i // 2 + 1] += self.de_f[i]
            print("雅可比矩阵")
            print(self.Jaccob)
            print("节点功率不平衡量")
            print(self.Y)
            print("节点电压修正量")
            print(self.de_f)
            print("节点电压新值")
            print(self.U)
            print("")
            # for j in range(2, 6):
            #     print("U%d" % j)
            #     printcomplex(self.U[j - 1])

    def PowerCal(self):#计算线路功率
        n=self.n
        # 计算线路功率
        print('平衡节点功率：')
        print("S1")
        S1 = self.U[0] * sum([self.YB[1 - 1][j].conjugate() * self.U[j].conjugate() for j in range(5)])
        printcomplex(S1)
        self.Sbranch=dict()
        #edge = [[1, 2], [2, 3], [1, 3], [3, 4], [2, 4], [2, 5], [4, 5]]
        edge=self.edge
        for i in range(1, n+1):
            for j in range(1, n+1):
                if [i, j] in edge or [j, i] in edge:
                    Sbranch=self.U[i - 1] * (self.U[i - 1].conjugate() - self.U[j - 1].conjugate()) * self.y[i - 1][j - 1].conjugate()
                    if Sbranch.real>0:
                        self.Sbranch[str(i)+','+str(j)]=printcomplex2(Sbranch)
                    print("S {} {}".format(i, j))
                    printcomplex(Sbranch)
        print("derta S sigma:")
        print(S1 + sum(self.S))

        self.S[0] = S1
        print("输电效率：")
        print("%.3f%%" % (sum([abs(self.S[i].real) for i in range(n) if i+1 not in self.inputNode]) / sum([abs(self.S[i].real) for i in range(n) if i+1 in self.inputNode]) * 100))

    def ultraV(self):
        print("最终电压")
        print(self.U)
        for u in self.U:
            rec2polar(u)

    def result(self):
        U=self.U
        S=self.S
        nodetype=[1 if i+1 in self.inputNode else 2 for i in range(self.n)]
        matrix=[[i+1,'%.4f'%abs(U[i]),'%.4f'%np.angle(U[i]),'%.4f'%S[i].real,'%.4f'%S[i].imag,nodetype[i]] for i in range(self.n)]
        df=pd.DataFrame(matrix,columns=['母线名','电压幅值','电压相角/(°)','有功功率','无功功率','类型(1发电机/2负荷)'])
        print(df)
        df.to_csv('Tide.txt',index=False)

    def G(self,i, j):
        return Re(self.YB[i - 1][j - 1])

    def B(self,i, j):
        return Im(self.YB[i - 1][j - 1])

    def f(self,i):
        return Im(self.U[i - 1])

    def e(self,i):
        return Re(self.U[i - 1])

    def P(self,i):
        ans = 0
        for j in range(1, self.n + 1):
            ans += self.e(i) * (self.G(i, j) * self.e(j) - self.B(i, j) * self.f(j)) + self.f(i) * (self.G(i, j) * self.f(j) + self.B(i, j) * self.e(j))
        return ans

    def Q(self,i):
        ans = 0
        for j in range(1, self.n + 1):
            ans += self.f(i) * (self.G(i, j) * self.e(j) - self.B(i, j) * self.f(j)) - self.e(i) * (self.G(i, j) * self.f(j) + self.B(i, j) * self.e(j))
        return ans
def norm(X):#范数
    n=len(X)
    ans=0
    for i in range(n):
        ans+=X[i]**2
    return ans**0.5
def Re(x):#实部
    return x.real
def Im(x):#虚部
    return x.imag
def printcomplex(x):#打印一定精度的复数
    if x.imag>=0:
        print("%.5f"%x.real,'+',"j%.5f"%x.imag)
    elif x.imag<0:
        print("%.5f"%x.real,'-',"j%.5f"%(-x.imag))
def printcomplex2(x):#返回一定精度的复数
    if x.imag>=0:
        return "%.4f"%x.real+'+'+"j%.4f"%x.imag
    elif x.imag<0:
        return "%.4f"%x.real+'-'+"j%.4f"%(-x.imag)
def rec2polar(x):#打印极坐标形式的一定精度的复数
    if x.real!=0:
        print("{:.4f}∠{:.4f}°".format(float(abs(x)),atan(x.imag/x.real)*180/pi))
    else:
        print("{:.4f}∠90°".format(float(abs(x))))
def rec2polar2(x):#返回极坐标形式的一定精度的复数
    if x.real!=0:
        return "{:.4f}∠{:.4f}°".format(float(abs(x)),atan(x.imag/x.real)*180/pi)
    else:
        return "{:.4f}∠90°".format(float(abs(x)))


