# This is a demo for LSM
# CaoYaqin, 2023/11/13
# caoyq23@mails.jlu.edu.cn

# Before you use LSM program, you should define the velocity file and head file
# This demo gives a program to generate the head file (program: makehead)
# Or you can use the data in data folder

# head file parameters
ntr=737
dtr=12.5
ftr=0
dftr=0

nshot=20
dshot=450
fshot=0

# forward parameters
naper=150
nt=2000
nx=737
nz=750
dt=0.004
dx=12.5
dz=8
amax=50
vfile=../../data/vel_marmousi # input 速度文件
reftfile=reflectivity.dat # output 初次计算的反射系数
orgfile=shot_0.dat # output 初次计算的地震数据 或用于对比的地震记录

# LSM parameters
R1name=resMig_0.dat
iters=60

../../bin/generate nx=$nx nz=$nz name=$R1name

# generate the head file 生成地震正演的观测系统文件
../../bin/makehead ntr=$ntr dtr=$dtr ftr=$ftr dftr=$dftr \
    nshot=$nshot dshot=$dshot fshot=$fshot 

# generate the seismic file 生成原始地震数据
../../bin/forwardGBM naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax \
    reftfile=$reftfile orgfile=$orgfile vfile=$vfile useOuterreflectivity=0 ntr=$ntr \
    dtr=$dtr nshot=$nshot

# copy the seismic file 复制一份原始地震数据用于对比
cp $orgfile orgseis.dat
orgseis=orgseis.dat

# LSM
for ((i=1;i<=$iters;i++))
do
    echo "\033[32mLSM iter: $i\033[0m"

    R=mig_$((i-1)).dat
    seis=shot_$((i-1)).dat
    newseis=shot_${i}.dat

    resMig=resMig_$((i-1)).dat
    newresMig=resMig_${i}.dat

    #forward 利用上一次的反射系数进行正演 如果是第一次迭代 就用原始的反射系数(全0)
    ../../bin/forwardGBM naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax \
            reftfile=$resMig orgfile=$seis vfile=$vfile useOuterreflectivity=1  ntr=$ntr \
            dtr=$dtr nshot=$nshot
        
    # subtract 用原始地震数据减去正演数据得到残差(梯度)
    ../../bin/subtract nx=$nx nz=$nz nt=$nt ntr=$ntr nshot=$nshot originShotname=$orgseis \
            newShotname=$seis subtractDataname=$newseis 

    # migration 用梯度进行偏移得到反射系数
    mpirun -n 4 ../../bin/migrationGBM naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz \
            amax=$amax seifile=$newseis vfile=$vfile imgfile=$R ntr=$ntr dtr=$dtr nshot=$nshot 

    # forward 利用新的反射系数进行正演
    ../../bin/forwardGBM naper=$naper nt=$nt nx=$nx nz=$nz dt=$dt dx=$dx dz=$dz amax=$amax \
            reftfile=$R orgfile=$seis vfile=$vfile useOuterreflectivity=1 ntr=$ntr dtr=$dtr \
            nshot=$nshot
    
    # LSM 最小二乘迭代函数 计算迭代步长 并更新反射系数
    ../../bin/lsm seisname=$seis migname=$R resMigname=$resMig newresMIGname=$newresMig \
            ntr=$ntr nshot=$nshot nx=$nx nz=$nz 
done
