# README

## 这是什么

计划是实现一个用于通信的rx dsp库，其中包含数据产生、信号同步提取、均衡算法和判决与误码计算等代码。

## 尚未完成

* 该readme的编写
* 新的均衡算法的实现
  * Volterra RLS FFE
  * Volterra 系数简化算法
  * MLSE（可能暂时基于MATLAB工具箱）
  * Linear DFE
  * Volterra DFE
  * 一些Machine Learning算法
* NRZ判决与误码计算
* 用Python实现

## 核心函数

###信号处理辅助函数 

```OriginalData = generateData(OriginalDataLength, PAM4Flag, NewPRBSGenerationFlag, SyncZerosLength)```

产生数据函数。该函数可产生源数据，供其他函数以使用。此外，该函数还可以产生用于载入PPG发送的数据。目前可选产生NRZ或者PAM4数据。

* 输入参数

  1. OriginalDataLength（可选）

     所产生的原始数据长度，默认为4096；

  2. PAM4Flag（可选）

     是否产生PAM4的标志：0为产生NRZ，1为产生PAM4。默认为1。

  3. NewPRBSGenerationFlag（可选）

     是否重新生成数据的标志：0为重新生成，1为读取./Original Data/Original_Data_4096.txt中的数据。默认为0。

  4. SyncZerosLength（可选）

     所产生给PPG的数据中同步部分0的个数，默认为50。

* 输出参数

  OriginalData 

  所产生或从文件中读取的NRZ或者PAM4源数据。

```[ExtractedSignalUS, OriginalSignalUS] = syncAndExtractSignal(SampledData, OriginalData, OverSamplingRatio, UpSamplingRatio)``` 

同步与信号提取函数。该函数可将所接收到的信号从DSO采集到的数据SampledData中同步并提取出来。

* 输入参数

  1. SampledData

     从DSO来的采样数据

  2. OriginalData

     原始发送数据

  3. OverSamplingRatio

     过采样率，通常等于DSO采样率/信号速率

  4. UpSamplingRatio（可选）

     输出信号上采样率。默认值为1，即为不进行上采样。

* 输出参数

  1. ExtractedSignalUS

     从SampledData中同步提取出来的信号，经过UpSamplingRatio倍上采样，是一个length(OriginalData)*UpSamplingRatio, 1的列向量。

  2. OriginalSignalUS

     经过UpSamplingRatio倍上采样的原始数据，也是一个length(OriginalData)*UpSamplingRatio, 1的列向量。

```[BitErrorRate, SymErrorRate, BitErrorNum] = decisionAndCalcBerPAM4(InputSignal, OriginalData, threshold)```

PAM4判决与误码计算函数。

* 输入参数

  1. InputSignal

     待判决的输入信号

  2. OriginalData

     用于计算误码的原始数据

  3. threshold（可选）

     PAM4判决电平，默认为[0.25; 0.5; 0.75]。

* 输出参数

  1. BitErrorRate

     误比特率。在计算误码率时，会计算错误一个符号时反转的比特数。

  2. SymErrorRate

     误符号率

  3. BitErrorNum

     误比特数

###均衡函数

```[output, w, costs] = linearFFEqualize(InputSignal, TrainingSignal, AlgType, FFETaps, alpha, epoch)```

linear FFE均衡函数。支持LMS和RLS两种算法。输入信号和训练信号会先被归一化到0-1。在目前的实现中所有的InputSignal都会被用于训练epoch次，在训练后会使用训练好的权重进行均衡。

* 输入参数

  1. InputSignal

     将被均衡的输入信号。

  2. TrainingSignal

     用于训练的目标序列。

  3. AlgType（可选）

     ‘lms’代表LMS算法，‘rls’代表RLS算法。默认为使用LMS算法。

  4. FFETaps（可选）

     FFE抽头数量，需为奇数，默认为5。

  5. alpha（可选）

     学习速率。默认为0.01。

  6. epoch（可选）

     训练次数。默认为1。

* 输出参数

  1. output

     均衡后的信号。一个length(InputSignal), 1的列向量。

  2. w

     训练好的权重。一个FFETaps, 1的列向量。

  3. costs

     每个训练epoch的cost向量，用于画收敛曲线以表现学习效果。


```[output, w, costs] = volLMSFFEqualize(InputSignal, TrainingSignal, chanLen, alpha1st, epoch, en2ndOrder, en3rdOrder, alpha2nd, alpha3rd)```

volterra FFE均衡函数。目前仅支持LMS算法，之后会增加对RLS算法的支持。输入信号和训练信号会先被归一化到0-1。在目前的实现中所有的InputSignal都会被用于训练epoch次，在训练后会使用训练好的权重进行均衡。

* 输入参数

  1. InputSignal

     将被均衡的输入信号。

  2. TrainingSignal

     用于训练的目标序列。

  3. chanLen（可选）

     信道长度。必须为奇数。默认为5。

  4. alpha1st（可选）

     volterra 1阶核的学习速率。默认为0.01。

  5. epoch（可选）

     训练次数。默认为1。

  6. en2ndOrder（可选）

     使用volterra 2阶核。默认为true。

  7. en3rdOrder（可选）

     使用volterra 3阶核。默认为true。

  8. alpha2nd（可选）

     volterra 2阶核的学习速率。默认与alpha1st相等。

  9. alpha3rd（可选）

     volterra 3阶核的学习速率。默认与alpha1st相等。

* 输出参数

  1. output

     均衡后的信号。一个length(InputSignal), 1的列向量。

  2. w

     训练好的权重。一个chanLen + kernel2ndSize + kernel3rdSize, 1的列向量。

  3. costs

     每个训练epoch的cost向量，用于画收敛曲线以表现学习效果。


## 一些例子

```LMSFFEqualizationExample.m```

这个例子展示了如何使用本库进行信号、同步、线性LMS均衡和判决的基本workflow。

```LMSFFEqualizationExample_2.m```

这个例子展示了如何对整个文件夹中所采集到的dat文件进行遍历以画出BER曲线。同样使用线性LMS均衡器。