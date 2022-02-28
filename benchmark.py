
import datetime
get=datetime.datetime.now

import cmath
import matplotlib.pyplot as plt
import time

#再現したい反応がないから改善するモチベ沸かない
class Polymer_FFT:
    def __init__(self,mol_length,volume,temperture):
        #FFT用
        self.pi = cmath.pi
        #反応記述用
        self.N=mol_length
        self.M=[0]*mol_length
        self.R=[0]*mol_length
        self.P=[0]*mol_length
        self.init=0

        self.T=temperture

        self.ki=self.Arrhenius(2.23*10**19/60,155.8*10**3).calc

    def drop(self,monomer,init):#途中でモノマーラジカルを入れたりすることを考えるとinitで入れるよりこうした方いいかもしれない 開始剤とモノマーで関数分けたら？
        self.M[1]+=monomer
        self.init+=init
        self.init_log=init

    def react(self,heatflow):
        self.T+=0.01
        print("温度:%f" % (self.T))
        print("開始剤:%f" % (self.init))
        delta_init=self.ki(self.T)*self.init
        self.init-=delta_init
        self.init=max(self.init,0)
        self.R[0]+=delta_init
        
        #反応する分子数を決める
        #1step分の時間幅がどれくらいにできるかはkiとkd次第かな 1stepで残ってるモノマー以上に反応されたらやばい
        #kdとktと連鎖移動からすべてのラジカルの行き先を決められるはず
        #将来的には開始剤ごとにラジカルを発生・モノマーごとに反応量を決める
        SUM_M,SUM_R=sum(self.M),sum(self.R)
        x=1/2*(SUM_R/(SUM_M+SUM_R))#アタックを受けて殺されるラジカルの割合
        R_die=[x*self.R[i] for i in range(self.N)]
        self.R=[self.R[i]-R_die[i] for i in range(self.N)]
        
        SUM_R_die=sum(R_die)
        M_per=[self.M[i]/(SUM_M) for i in range(self.N)]
        R_die_per=[R_die[i]/(SUM_R_die) for i in range(self.N)]
        print("R_die_per:%f" % (sum(R_die_per)))

        RtoM=[self.R[i]-R_die[i] for i in range(self.N)]
        RtoR=R_die[:]

        #FFTする
        delta_P=(self.convolution(R_die_per,RtoR))
        delta_R=(self.convolution(M_per,RtoM))
        #delta_R=[0]+RtoM[:-1]

        """
        dead_N=len(delta_P)-len(C)
        dead+=sum([delta_C[i]*i for i in range(N,N+dead_N)])+sum([delta_B[i]*i for i in range(N,N+dead_N)])
        """
        print("残りモノマー:%f ラジカル反応数:%f 死ぬラジカル数%f 生まれたポリマー数%f" % (SUM_M,SUM_R,SUM_R_die,sum(delta_P)))#死ぬラジカル数とポリマー数が一緒じゃないとおかしくない？
        #print("計算上死んだモノマー:%f" % (dead))
        
        #モノマー
        self.M[1]-=(SUM_R-2*SUM_R_die)*self.M[1]/sum(self.M)
        #ラジカル
        self.R=[delta_R[i] for i in range(self.N)]
        #ポリマー
        self.P=[self.P[i]+delta_P[i] for i in range(self.N)]
        
        print('ポリマー量:%f' % (sum([self.M[i]*i for i in range(self.N)])+sum([self.R[i]*i for i in range(self.N)])+sum([self.P[i]*i for i in range(self.N)])))

        print('残りモノマー数:%f' % (sum(self.M)))


        return
    
    def display(self):
        fig=plt.figure()
        ax1=fig.add_subplot(1, 1, 1)
        ax1.plot([self.P[i]*i for i in range(self.N)])
        ax1.plot([self.R[i]*i for i in range(self.N)],color="red")
        ax1.set_ylim(0,3000)
        ax1.set_xlim(0,self.N)
        plt.title("ki=%f,init=%i" %(self.ki(self.T),self.init_log))
        plt.show()

        return
    
    class Arrhenius:
        def __init__(self,A,Ea):
            self.A=A
            self.Ea=Ea
            self.R=8.314
        
        def calc(self,T):#ki(T=180)とかできるように
            #return self.A*cmath.exp(-self.Ea/(self.R*(T+273.15))).real
            return 0.1
    def convolution(self, g, h):
        l = len(g) + len(h)
        n = 1 << l.bit_length()
        
        g = g + [0] * (n - len(g))
        h = h + [0] * (n - len(h))
       
        gg = self.fft(g, n)
        hh = self.fft(h, n)
        ff = [gg[i] * hh[i] for i in range(n)]
 
        return self.inv_fft(ff, n)
 
    def fft(self, f, n, sgn = 1):
        if n == 1: return f
 
        f0 = [f[2 * i] for i in range(n // 2)]
        f1 = [f[2 * i + 1] for i in range(n // 2)]
        f0 = self.fft(f0, n // 2, sgn)
        f1 = self.fft(f1, n // 2, sgn)
 
        zeta = cmath.rect(1, sgn * self.pi * 2 / n)
        pow_zeta = 1
 
        for i in range(n):
            f[i] = f0[i % (n // 2)] + pow_zeta * f1[i % (n // 2)]
            pow_zeta *= zeta
 
        return f
 
    def inv_fft(self, f, n):
        ret = self.fft(f, n, -1)
        for i in range(n):
            ret[i] = (ret[i].real / n)
 
        return ret

test=Polymer_FFT(5*10**3,2,160)
test.drop(10**6,10*10**3)

count=0
start=get()
while True:
    print(get()-start)
    count+=1
    test.react(-1)
    if test.M[1]<10**6/100:
        test.display()
        print("----終了----")
        print(count)
        break

import cmath
import matplotlib.pyplot as plt
import time

#再現したい反応がないから改善するモチベ沸かない
class Polymer_FFT:
    def __init__(self,mol_length,volume,temperture):
        #FFT用
        self.pi = cmath.pi
        #反応記述用
        self.N=mol_length
        self.M=[0]*mol_length
        self.R=[0]*mol_length
        self.P=[0]*mol_length
        self.init=0

        self.T=temperture

        self.ki=self.Arrhenius(2.23*10**19/60,155.8*10**3).calc

    def drop(self,monomer,init):#途中でモノマーラジカルを入れたりすることを考えるとinitで入れるよりこうした方いいかもしれない 開始剤とモノマーで関数分けたら？
        self.M[1]+=monomer
        self.init+=init
        self.init_log=init

    def react(self,heatflow):
        self.T+=0.01
        print("温度:%f" % (self.T))
        print("開始剤:%f" % (self.init))
        delta_init=self.ki(self.T)*self.init
        self.init-=delta_init
        self.init=max(self.init,0)
        self.R[0]+=delta_init
        
        #反応する分子数を決める
        #1step分の時間幅がどれくらいにできるかはkiとkd次第かな 1stepで残ってるモノマー以上に反応されたらやばい
        #kdとktと連鎖移動からすべてのラジカルの行き先を決められるはず
        #将来的には開始剤ごとにラジカルを発生・モノマーごとに反応量を決める
        SUM_M,SUM_R=sum(self.M),sum(self.R)
        x=1/2*(SUM_R/(SUM_M+SUM_R))#アタックを受けて殺されるラジカルの割合
        R_die=[x*self.R[i] for i in range(self.N)]
        self.R=[self.R[i]-R_die[i] for i in range(self.N)]
        
        SUM_R_die=sum(R_die)
        M_per=[self.M[i]/(SUM_M) for i in range(self.N)]
        R_die_per=[R_die[i]/(SUM_R_die) for i in range(self.N)]
        print("R_die_per:%f" % (sum(R_die_per)))

        RtoM=[self.R[i]-R_die[i] for i in range(self.N)]
        RtoR=R_die[:]

        #FFTする
        delta_P=(self.convolution(R_die_per,RtoR))
        delta_R=(self.convolution(M_per,RtoM))
        
        #delta_R=[0]+RtoM[:-1]

        """
        dead_N=len(delta_P)-len(C)
        dead+=sum([delta_C[i]*i for i in range(N,N+dead_N)])+sum([delta_B[i]*i for i in range(N,N+dead_N)])
        """
        print("残りモノマー:%f ラジカル反応数:%f 死ぬラジカル数%f 生まれたポリマー数%f" % (SUM_M,SUM_R,SUM_R_die,sum(delta_P)))#死ぬラジカル数とポリマー数が一緒じゃないとおかしくない？
        #print("計算上死んだモノマー:%f" % (dead))
        
        #モノマー
        self.M[1]-=(SUM_R-2*SUM_R_die)*self.M[1]/sum(self.M)
        #ラジカル
        self.R=[delta_R[i] for i in range(self.N)]
        #ポリマー
        self.P=[self.P[i]+delta_P[i] for i in range(self.N)]
        
        print('ポリマー量:%f' % (sum([self.M[i]*i for i in range(self.N)])+sum([self.R[i]*i for i in range(self.N)])+sum([self.P[i]*i for i in range(self.N)])))

        print('残りモノマー数:%f' % (sum(self.M)))


        return
    
    def display(self):
        fig=plt.figure()
        ax1=fig.add_subplot(1, 1, 1)
        ax1.plot([self.P[i]*i for i in range(self.N)])
        ax1.plot([self.R[i]*i for i in range(self.N)],color="red")
        ax1.set_ylim(0,3000)
        ax1.set_xlim(0,self.N)
        plt.title("ki=%f,init=%i" %(self.ki(self.T),self.init_log))
        plt.show()

        return
    
    class Arrhenius:
        def __init__(self,A,Ea):
            self.A=A
            self.Ea=Ea
            self.R=8.314
        
        def calc(self,T):#ki(T=180)とかできるように
            #return self.A*cmath.exp(-self.Ea/(self.R*(T+273.15))).real
            return 0.1
    def convolution(self, g, h):
        l = len(g) + len(h)
        n = 1 << l.bit_length()
        
        g = g + [0] * (n - len(g))
        h = h + [0] * (n - len(h))
       
        gg = self.fft(g, n)
        hh = self.fft(h, n)
        ff = [gg[i] * hh[i] for i in range(n)]
 
        return self.inv_fft(ff, n)
 
    def fft(self, f, n, sgn = 1):
        if n == 1: return f
 
        f0 = [f[2 * i] for i in range(n // 2)]
        f1 = [f[2 * i + 1] for i in range(n // 2)]
        f0 = self.fft(f0, n // 2, sgn)
        f1 = self.fft(f1, n // 2, sgn)
 
        zeta = cmath.rect(1, sgn * self.pi * 2 / n)
        pow_zeta = 1
 
        for i in range(n):
            f[i] = f0[i % (n // 2)] + pow_zeta * f1[i % (n // 2)]
            pow_zeta *= zeta
 
        return f
 
    def inv_fft(self, f, n):
        ret = self.fft(f, n, -1)
        for i in range(n):
            ret[i] = (ret[i].real / n)
 
        return ret

test=Polymer_FFT(3*10**4,2,160)
test.drop(10**6,10*10**3)

count=0
start=get()
while True:
    print(get()-start)
    count+=1
    test.react(-1)
    if test.M[1]<10**6/100:
        test.display()
        print("----終了----")
        print(count)
        break

import cmath
import matplotlib.pyplot as plt
import time

#再現したい反応がないから改善するモチベ沸かない
class Polymer_FFT:
    def __init__(self,mol_length,volume,temperture):
        #FFT用
        self.pi = cmath.pi
        #反応記述用
        self.N=mol_length
        self.M=[0]*mol_length
        self.R=[0]*mol_length
        self.P=[0]*mol_length
        self.init=0

        self.T=temperture

        self.ki=self.Arrhenius(2.23*10**19/60,155.8*10**3).calc

    def drop(self,monomer,init):#途中でモノマーラジカルを入れたりすることを考えるとinitで入れるよりこうした方いいかもしれない 開始剤とモノマーで関数分けたら？
        self.M[1]+=monomer
        self.init+=init
        self.init_log=init

    def react(self,heatflow):
        self.T+=0.01
        print("温度:%f" % (self.T))
        print("開始剤:%f" % (self.init))
        delta_init=self.ki(self.T)*self.init
        self.init-=delta_init
        self.init=max(self.init,0)
        self.R[0]+=delta_init
        
        #反応する分子数を決める
        #1step分の時間幅がどれくらいにできるかはkiとkd次第かな 1stepで残ってるモノマー以上に反応されたらやばい
        #kdとktと連鎖移動からすべてのラジカルの行き先を決められるはず
        #将来的には開始剤ごとにラジカルを発生・モノマーごとに反応量を決める
        SUM_M,SUM_R=sum(self.M),sum(self.R)
        x=1/2*(SUM_R/(SUM_M+SUM_R))#アタックを受けて殺されるラジカルの割合
        R_die=[x*self.R[i] for i in range(self.N)]
        self.R=[self.R[i]-R_die[i] for i in range(self.N)]
        
        SUM_R_die=sum(R_die)
        M_per=[self.M[i]/(SUM_M) for i in range(self.N)]
        R_die_per=[R_die[i]/(SUM_R_die) for i in range(self.N)]
        print("R_die_per:%f" % (sum(R_die_per)))

        RtoM=[self.R[i]-R_die[i] for i in range(self.N)]
        RtoR=R_die[:]

        """
        #FFTする
        delta_P=(self.convolution(R_die_per,RtoR))
        delta_R=(self.convolution(M_per,RtoM))
        """
        delta_P=[0 for i in range(self.N)]
        for i in range(len(R_die_per)-1):
            for j in range(i,len(RtoR)-i):
                delta_P[i+j]+=R_die_per[i]*RtoR[j]
        delta_R=[0 for i in range(self.N)]
        for i in range(len(M_per)-1):
            for j in range(i,len(RtoM)-i):
                delta_R[i+j]+=M_per[i]*RtoM[j]

        #delta_R=[0]+RtoM[:-1]

        """
        dead_N=len(delta_P)-len(C)
        dead+=sum([delta_C[i]*i for i in range(N,N+dead_N)])+sum([delta_B[i]*i for i in range(N,N+dead_N)])
        """
        print("残りモノマー:%f ラジカル反応数:%f 死ぬラジカル数%f 生まれたポリマー数%f" % (SUM_M,SUM_R,SUM_R_die,sum(delta_P)))#死ぬラジカル数とポリマー数が一緒じゃないとおかしくない？
        #print("計算上死んだモノマー:%f" % (dead))
        
        #モノマー
        self.M[1]-=(SUM_R-2*SUM_R_die)*self.M[1]/sum(self.M)
        #ラジカル
        self.R=[delta_R[i] for i in range(self.N)]
        #ポリマー
        self.P=[self.P[i]+delta_P[i] for i in range(self.N)]
        
        print('ポリマー量:%f' % (sum([self.M[i]*i for i in range(self.N)])+sum([self.R[i]*i for i in range(self.N)])+sum([self.P[i]*i for i in range(self.N)])))

        print('残りモノマー数:%f' % (sum(self.M)))


        return
    
    def display(self):
        fig=plt.figure()
        ax1=fig.add_subplot(1, 1, 1)
        ax1.plot([self.P[i]*i for i in range(self.N)])
        ax1.plot([self.R[i]*i for i in range(self.N)],color="red")
        ax1.set_ylim(0,3000)
        ax1.set_xlim(0,self.N)
        plt.title("ki=%f,init=%i" %(self.ki(self.T),self.init_log))
        plt.show()

        return
    
    class Arrhenius:
        def __init__(self,A,Ea):
            self.A=A
            self.Ea=Ea
            self.R=8.314
        
        def calc(self,T):#ki(T=180)とかできるように
            #return self.A*cmath.exp(-self.Ea/(self.R*(T+273.15))).real
            return 0.1
    def convolution(self, g, h):
        l = len(g) + len(h)
        n = 1 << l.bit_length()
        
        g = g + [0] * (n - len(g))
        h = h + [0] * (n - len(h))
       
        gg = self.fft(g, n)
        hh = self.fft(h, n)
        ff = [gg[i] * hh[i] for i in range(n)]
 
        return self.inv_fft(ff, n)
 
    def fft(self, f, n, sgn = 1):
        if n == 1: return f
 
        f0 = [f[2 * i] for i in range(n // 2)]
        f1 = [f[2 * i + 1] for i in range(n // 2)]
        f0 = self.fft(f0, n // 2, sgn)
        f1 = self.fft(f1, n // 2, sgn)
 
        zeta = cmath.rect(1, sgn * self.pi * 2 / n)
        pow_zeta = 1
 
        for i in range(n):
            f[i] = f0[i % (n // 2)] + pow_zeta * f1[i % (n // 2)]
            pow_zeta *= zeta
 
        return f
 
    def inv_fft(self, f, n):
        ret = self.fft(f, n, -1)
        for i in range(n):
            ret[i] = (ret[i].real / n)
 
        return ret

test=Polymer_FFT(3*10**4,2,160)
test.drop(10**6,10*10**3)

count=0
start=get()
while True:
    print(get()-start)
    count+=1
    test.react(-1)
    test.display()
    if test.M[1]<10**6/100:
        test.display()
        print("----終了----")
        print(count)
        break