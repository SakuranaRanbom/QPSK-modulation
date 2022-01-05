from math import log
from numpy import random,fft
from scipy.fftpack import fft,ifft,fftshift
from numpy.core.fromnumeric import size
from numpy.fft.helper import fftshift

from scipy import signal, special
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker 
import gausnoise
def wgn(x, snr):
    # snr = 10**(snr/10.0)
    # xpower = np.sum(x**2)/len(x)
    # npower = xpower / snr
    # return np.random.randn(len(x)) * np.sqrt(npower)
    Ps = np.sum(abs(x)**2)/len(x)
    Pn = Ps/(10**(snr/10))
    return np.random.randn(len(x)) * np.sqrt(Pn)

def bitxor(lst1, lst2):
    """两个二进制列表的逐位异或"""
    if len(lst1) != len(lst2):
        print('两个列表长度不一致！')
        return []
    else:
        return [lst1[x]^lst2[x] for x in range(len(lst1))]

def cm_sm33(snr_in_dB):
    N=100
    E=1

    snr=pow(10,snr_in_dB/10)
    sgma=np.sqrt(E/snr)/2
    s00=[1,0]
    s01=[0,1]
    s11=[-1,0]
    s10=[0,-1]
    desource1 = np.zeros(N)
    desource2 = np.zeros(N)
    numofsymbolerror = 0
    numofbiterror=0
    counter=0
    while numofbiterror<100:
        for i in range(1,N):
            temp=np.random.rand()
            if temp<0.25:
                desource1[i]=0
                desource2[i]=0
            elif temp<0.5:
                desource1[i]=0
                desource2[i]=1
            elif temp<0.75:
                desource1[i]=1
                desource2[i]=0
            else:
                desource1[i]=1
                desource2[i]=1
        for i in range(1,N):
            n=sgma*random.randn(1,2)
            if (desource1[i]==0) and (desource2[i]==0):
                r=s00+n
            elif (desource1[i]==0) and (desource2[i]==1):
                r=s01+n
            elif (desource1[i]==1) and (desource2[i]==0):
                r=s10+n
            else:
                r=s11+n
            c00=np.dot(r,s00)
            c01=np.dot(r,s01)
            c10=np.dot(r,s10)
            c11=np.dot(r,s11)
            c_max=max([c00,c01,c11,c10])
            if c00==c_max:
                decis1=0
                decis2=0
            elif c01==c_max:
                decis1=0
                decis2=1
            elif c11==c_max:
                decis1=1
                decis2=1
            else:
                decis1=1
                decis2=0
            symbolerror=0
            if decis1!=desource1[i]:
                numofbiterror=numofbiterror+1
                symbolerror=1
            if decis2!=desource2[i]:
                numofbiterror=numofbiterror+1
                symbolerror=1
            if symbolerror==1:
                numofsymbolerror=numofsymbolerror+1
    
        counter=counter+1
    return numofbiterror






T = 1               #基带信号宽度，也就是频率
nb = 2000       #定义传输的比特数
delta_T = T/200     #采样间隔
fs = 1/delta_T      #采样频率
fc = 10/T           #载波频率
SNR = [-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10]             #信噪比
errorcnt = []
t = np.arange(0, nb*T, delta_T)
N = len(t)

# 产生基带信号
data = [1 if x > 0.5 else 0 for x in np.random.randn(1, nb)[0]]  #调用随机函数产生任意在0到1的1*nb的矩阵，大于0.5显示为1，小于0.5显示为0
data0 = []                             #创建一个1*nb/delta_T的零矩阵
for q in range(nb):
    data0 += [data[q]]*int(1/delta_T)  #将基带信号变换成对应波形信号
# plt.figure(12,20)
# plt.plot(data)
# plt.show()
# 调制信号的产生
data1 = []      #创建一个1*nb/delta_T的零矩阵
datanrz = np.array(data)*2-1              #将基带信号转换成极性码,映射
for q in range(nb):
    data1 += [datanrz[q]]*int(1/delta_T)  #将极性码变成对应的波形信号
    
idata = datanrz[0:(nb-1):2]       #串并转换，将奇偶位分开，间隔为2，i是奇位 q是偶位
qdata = datanrz[1:nb:2]         
ich = []                          #创建一个1*nb/delta_T/2的零矩阵，以便后面存放奇偶位数据
qch = []         
for i in range(int(nb/2)):
    ich += [idata[i]]*int(1/delta_T)    #奇位码元转换为对应的波形信号
    qch += [qdata[i]]*int(1/delta_T)    #偶位码元转换为对应的波形信号

a = []     #余弦函数载波
b = []     #正弦函数载波
for j in range(int(N/2)):
    a.append(np.math.sqrt(2/T)*np.math.cos(2*np.math.pi*fc*t[j]))    #余弦函数载波
    b.append(np.math.sqrt(2/T)*np.math.sin(2*np.math.pi*fc*t[j]))    #正弦函数载波
idata1 = np.array(ich)*np.array(a)          #奇数位数据与余弦函数相乘，得到一路的调制信号
qdata1 = np.array(qch)*np.array(b)          #偶数位数据与余弦函数相乘，得到另一路的调制信号
s = idata1 + qdata1      #将奇偶位数据合并，s即为QPSK调制信号

# plt.figure(figsize=(30,40))
# plt.subplot(1,1,1)
# plt.plot(datanrz)
# plt.title('origin')
# plt.figure(figsize=(30,40))
# plt.subplot(1,1,1)
# plt.plot(idata)
# plt.title('I')
# plt.figure(figsize=(30,40))
# plt.subplot(1,1,1)
# plt.plot(qdata)
# plt.title('Q')


# plt.figure(figsize=(14,12))
# plt.subplot(3,1,1)
# plt.plot(idata1)
# plt.title('I')
# plt.axis([0,500,-3,3])
# plt.subplot(3,1,2)
# plt.plot(qdata1)
# plt.title('Q')
# plt.axis([0,500,-3,3])
# plt.subplot(3,1,3)
# plt.plot(s)
# plt.title('QPSK')
# plt.axis([0,500,-3,3])
# plt.show()
# plt.figure(figsize=(30,40))
# s_f=abs(fftshift(fft(s)))
# plt.plot(s_f)
# plt.title('QPSK-f')
# plt.show()

#调制信号经过高斯信道
#s11 = wgn(s, SNR)     #高斯噪声曲线
#s11 = 
#noise =np.random.rand(s.shape[0], s)

#s1 = s + s11          #加上高斯信道之后的信号 
for ij in range(len(SNR)):
    ta = a
    tb =b
    s1 = gausnoise.AWGN(s,SNR[ij],1)
    print(len(s1))
    print(len(a))
    #解调（高斯信道）
    idata2 = s1*np.array(ta)       #这里面其实隐藏了一个串并转换的过程
    qdata2 = s1*np.array(tb)       #对应的信号与正余弦信号相乘
    demodata0 = []
    [tb, ta] = signal.butter(2, 2*fc/fs)
    idata22 = signal.filtfilt(tb, ta, idata2)
    qdata22 = signal.filtfilt(tb, ta, qdata2)
    for i in range(1,size(idata22)):
        demodata0.append(idata22[i])
        demodata0.append(qdata22[i])
    #demodata0 = idata22 + qdata22

    #print(idata22)
    #print(qdata22)
    #print(demodata0)

    # plt.figure(figsize=(30,40))
    # plt.subplot(1,1,1)
    # plt.plot(s1)
    # plt.title('Awgn')
    # plt.figure(figsize=(30,40))
    # plt.plot(idata1)
    # plt.plot(idata2)
    # plt.figure(figsize=(30,40))

    # plt.subplot(1, 1, 1)
    # plt.plot(demodata0)

    # plt.title('dismodel output')
    # plt.axis([0, 5000, -4, 4])
    # plt.show()

    Awgn_ichsum = []
    Awgn_qchsum = []
    for i in range(int(nb/2)):
        Awgn_ichsum.append(np.sum(idata22[i*int(1/delta_T):(i+1)*int(1/delta_T)])*delta_T)
        Awgn_qchsum.append(np.sum(qdata22[i*int(1/delta_T):(i+1)*int(1/delta_T)])*delta_T)

    # plt.figure(figsize=(12,10))
    # for i in range(int(nb/2)):
    #     plt.scatter(idata[i], qdata[i], s=150, c='r', marker='+')  #散点图，s表示点型的大小，c表示颜色，marker表示点型
    #     plt.title('constellation diagram')
    #     plt.axis([-2, 2, -2, 2])
    #     plt.scatter(Awgn_ichsum[i], Awgn_qchsum[i], s=150, c='black', marker='*')
    #     plt.legend(['theoretical','actual'])
    # plt.show()


    #解调后从时域积分（求和）判决
    idata3 = []            #建立1*nb/2数组，以存放判决之后的奇信号
    qdata3 = []            #建立1*nb/2数组，以存放判决之后的偶信号
    #抽样判决的过程，与0作比较，data>=0,则置1，否则置0
    for i in range(int(nb/2)):
        if np.sum(idata22[i*int(1/delta_T):(i+1)*int(1/delta_T)]) >= 0:
            idata3.append(1)
        else:
            idata3.append(0)
        if np.sum(qdata22[i*int(1/delta_T):(i+1)*int(1/delta_T)]) >= 0:
            qdata3.append(1)
        else:
            qdata3.append(0)

    #将判决后的数据存放进数据组
    demodata = []
    for i in range(nb):
        if i % 2 == 0:
            demodata.append(idata3[i//2])   #并串变换，存放奇数位
        else:
            demodata.append(qdata3[i//2])   #并串变换，存放偶数位
            
    #为了显示，将它变成波形信号（即传输一个1代表单位宽度的高电平）
    demodata1 = []          #创建一个1*nb/delta_T的零矩阵
    for i in range(nb):
        demodata1 += [demodata[i]]*int(1/delta_T)    #将极性码变成对应的波形信号

    error1 = bitxor(data, demodata)
    error11 = np.sum(error1)
    Awgn_num_BER = error11/nb
    errorcnt.append(error11)



print(errorcnt)
#误码率曲线
SNRindB1 = np.arange(0, 7, 1)
SNRindB2 = np.arange(0, 7, 0.1)
#print(data)
#print(demodata)
#高斯信道
smld_bit_awgn_err_prb = []
for i in range(len(SNRindB1)):
    pb1 = errorcnt[i+2]/2000.0

    smld_bit_awgn_err_prb.append(pb1)
    
#理论曲线
theo_err_awgn_prb = []
for i in range(len(SNRindB2)):
    SNR = np.exp(SNRindB2[i]*np.log(10)/10)                        #信噪比
    theo_err_awgn_prb.append(0.5*special.erfc(np.math.sqrt(SNR)))  #高斯噪声理论误码率
plt.figure(figsize=(12, 8))
plt.semilogy(SNRindB2, theo_err_awgn_prb, 'b')
plt.semilogy(SNRindB1, smld_bit_awgn_err_prb, 'r*')
plt.title('Ber-Eb/No figure')
plt.axis([0, 6, 1.0e-3, 1e-1])
plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0e'))
plt.xlabel('Eb/No')
plt.ylabel('BER')
plt.legend(['Theoretical AWGN','Actual AWGN',])
plt.show()
print('error11:')
print(error11)