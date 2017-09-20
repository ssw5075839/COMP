# Compressed-Sensing Sigma Delta ADC with SPG, OMP and CoSaMP
This is a MATLAB simulation program for compressed-sensing sigma-delta analog to digital converter (ADC). For compressed sensing theory, please refer to my file [Theoretical_Background.pdf](https://github.com/ssw5075839/COMP/blob/master/Theoretical_Background.pdf). In this Readme file, I will only briefly introduce the concepts and then give the code running instruction.

## Compressed Sensing in Sigma Delta Analog to Digital Converter (ADC)
We have covered all the mathematical fundamentals to understand compressed sensing. We know that when A obeys RIP, compressed sensing is feasible and the reconstruction algorithm is available either through BP or BPDN problem. Now, let us take a look on how to integrate this process into Sigma Delta ADC.

The basic concepts of sigma delta ADC can be found in this book [Analog Integrated Circuit Design](https://www.amazon.com/Analog-Integrated-Circuit-Design-David/dp/0471144487/). Now let us analyze circuits model and explain how it works.

Below is an traditional sigma delta ADC:<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/normal_ds.PNG"></p>X(t) stands for input signal, Y(t) stands for ADC output, Res(t) stands for the residual value after the first subtractor, E(t) stands for the quantization noise introduced by the single bit quantizer and Int(t) stands for integrator output.

From this graph we can get three equations:<p align="center"><img src="http://latex.codecogs.com/gif.latex?\begin{cases}Y(t)=Int(t)+E(t)\\Int(t)=Res(t)+Int(t-1)\\Res(t)=X(t)-Y(t-1)\end{cases}" border="0" align="center"/></p>

It is not hard to solve it and get<p align="center"><img src="http://latex.codecogs.com/gif.latex?Y(t)=X(t)+E(t)-E(t-1)" border="0" align="center"/></p>

Here we get the behaviour of the traditional sigma delta ADC: since quantization noise will be deducted in the long rung (more details will be shown later on), Y(t) is approximately eqaul to input data X(t), namely, convert analog input signal X(t) to digital signal Y(t). This is achieved by adding only shapped quantization noise <img src="http://latex.codecogs.com/gif.latex?E(t)-E(t-1)" border="0" align="center"/>. 

Below is the compressed sensing sigma delta ADC structure used in this project:<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/cs_ds.PNG"></p>

The symbology is the same as the traditional one above except for two multipliers and <img src="http://latex.codecogs.com/gif.latex?\Phi" border="0" align="center"/>. As we discussed in the theoretical part, <img src="http://latex.codecogs.com/gif.latex?\Phi" border="0" align="center"/> is +/- 1 sampled from Bernoulli distribution. In the real circuits, this part will be a common random number generator. With the multipliers and <img src="http://latex.codecogs.com/gif.latex?\Phi" border="0" align="center"/>, the equations for sigma delta ADC changes to:<p align="center"><img src="http://latex.codecogs.com/gif.latex?\begin{cases}Y(t)=Int(t)+E(t)\\Int(t)=Res(t)+Int(t-1)\\Res(t)=\Phi(t)(X(t)-\Phi(t)Y(t-1))\end{cases}" border="0" align="center"/></p>

Please note that <img src="http://latex.codecogs.com/gif.latex?\Phi" border="0" align="center"/> is either +1 or -1 so <img src="http://latex.codecogs.com/gif.latex?\Phi^2" border="0" align="center"/> is always 1. Then, it is not hard to solve the equations and get<p align="center"><img src="http://latex.codecogs.com/gif.latex?Y(t)=\Phi(t)X(t)+E(t)-E(t-1)" border="0" align="center"/></p>

Compare to the traditional sigma delta ADC, the new one still has the noise shaping capability. But now it is doing compressed sensing. Let us consider the setup below. For each of the sigma delta ADC, we put an adder after Y(t) to accumulate Y(t), namely <img src="http://latex.codecogs.com/gif.latex?\mu=\sum_{t=1}^{n}Y(t)" border="0" align="center"/>, where n is the total sampling point of input signal X(t). Then we deploy m such ADCs in parallel:<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/mpaths.PNG" border="0" align="center"/></p>

For each path of ADC, use the equatio we just get we can derive that:<p align="center"><img src="http://latex.codecogs.com/svg.latex?\vec{\mu%20}=\begin{bmatrix}\mu%20_1\\...\\\mu_m%20\end{bmatrix}=\begin{bmatrix}\sum_{t=1}^{n}Y_1(t)\\...\\\sum_{t=1}^{n}Y_m(t)\end{bmatrix}=\begin{bmatrix}\sum_{t=1}^{n}\Phi_1(t)X(t)\\...\\\sum_{t=1}^{n}\Phi_m(t)X(t)\end{bmatrix}+\begin{bmatrix}E_1(n)-E_1(1)\\...\\E_m(n)-E_m(1)\end{bmatrix}\\\\\\=\Phi%20X+\sigma" border="0" align="center"/></p>

Here, we treat each path's accumulated output <img src="http://latex.codecogs.com/svg.latex?\mu_i" border="0" align="center"/> as one element of the output vector <img src="http://latex.codecogs.com/svg.latex?\vec{\mu}" border="0" align="center"/>. Each path has its own random number generator so <img src="http://latex.codecogs.com/svg.latex?\Phi=\begin{bmatrix}\Phi_1(1)%20&%20...%20&%20\Phi_1(n)\\...&%20&...\\\Phi_m(1)&...&\Phi_m(n)\end{bmatrix}" border="0" align="center"/> is a random matrix obeying RIP property. X(t) is the input signal applied to all the m paths and is sampled at time from 1 to n. <img src="http://latex.codecogs.com/svg.latex?\sigma=\begin{bmatrix}E_1(n)-E_1(1)\\...\\E_m(n)-E_m(1)\end{bmatrix}" border="0" align="center"/> is the accumulated shaped quantization noise. As stated before, shaping noise around DC is a important behavior of sigma delta ADC. When quantization noise is around DC, noise value should approach its average in the long run and therefore <img src="http://latex.codecogs.com/svg.latex?E_i(n)\approx%20E_i(1)" border="0" align="center"/>. With a low path digital filter, we can almost eliminate this shaped quantization noise.

By grouping the m paths accumulated outputs together, we get an interesting result:<p align="center"><img src="http://latex.codecogs.com/svg.latex?\vec{\mu}=\Phi%20X+\sigma" border="0" align="center"/></p>

Note that this is exactly the basis pursuit denoising (BPDN) problem as we stated at the beginning. With accumulated output <img src="http://latex.codecogs.com/svg.latex?\vec{\mu}" border="0" align="center"/> and random matrix <img src="http://latex.codecogs.com/svg.latex?\Phi" border="0" align="center"/>, we can use either SPGL1 or OMP or CoSaMP to reconstruct the sparse input signal X(t). Note that X(t) has a total length of n (sample points) and <img src="http://latex.codecogs.com/svg.latex?\vec{\mu}" border="0" align="center"/> has only m elements. Therefore, we compressed a sparse signal with length n to m. As stated in the theoretical part, as long as m is greater than the order of magnitude of <img src="http://latex.codecogs.com/svg.latex?S\log%20n" border="0" align="center"/>, the reconstruction is feasible. By transmitting the m symbols instead of n symbols, we could save a lot of wireless bandwidth resourses.

In the real world, most nature signal is not sparse in time domain. However, they are sparse in frequency domain. For example, a pure sine wave <img src="http://latex.codecogs.com/svg.latex?X(t)%20=%20sin(2\pi%20f\cdot%20t)" border="0" align="center"/> is infinitely long in time domain, but is only a spike at f in frequency domain. To really compressed sensing the nature signal, the time domain signal will be multiplied by an [FFT matrix](https://en.wikipedia.org/wiki/Fast_Fourier_transform) to convert to frequency domain. Then the frequency-sparse signal will be compressed and reconstructed. Finally, multiply the reconstruction with the inverse of FFT matrix to get the original signal in time domain. Sigma Delta ADC naturally did the FFT conversion and we only need to do the inverse FFT conversion in MATLAB.

## Code Running Instruction
In this matlab demo of compressed sensing sigma delta ADC, we will implment a one-bit first order sigma delta ADC behaviour model together with compressed sensing. Then we will feed ADC outputs in one of the SPGL1, OMP or CoSaMP to try to reconstruct sparse signal. If the signal is sparse enough in frequency domain and ADC path number is large enough, then all the algorithm should be able to reconstruct the original signal.

Please change direcotry to COMP/codes and start MATLAB from there (optionally you can add this directory to your MATLAB path). Run spgsetup.m and then we can start to run demo.

To run the code, first define a MATLAB structure opts. (This step is optional. If you don't define opts, the code will use all the default settings.) To define a struct, use command ```opts = struct('field1',value1,...,'fieldN',valueN)```. Each option comes with a ('field', value) pair. The available filed name and suggested value will be discussed in the next section. If you don't define some name field, the default setting of that option will be used.

Then run the function with command ```SD_COMP_DEMO(opts)``` (or simply ```SD_COMP_DEMO(struct())``` if you prefer to use default settings). The demo will output the reconstruction algorithm running status in the console and the oupt three graphs to compare the original signal agains reconstructed signal in detail, both in time domain and frequency domain.

## Explanation of Parameters in Program
The available option fields and its suggested values are:

1. 'PATHS'
<br />This is the setting for how many ADCs to use for compressed sensing. It is equal to the 'm' we mentioend above and in the theoretical background. Default is 15 which is suitable for recovering 2 to 3 sparse signals.

2. 'DATA_LENGTH'
<br />This is the setting for how many sampling points to use for the time domain signal. The larger it is, the more details the signal is recorded in. It is equal to the 'n' we mentioend above and in the theoretical background. Default is 256.

3. 'WARMUP_CYCLES'
<br />How many warm up cycles for ADC before it starts compressed sensing job. This is used to eleminate the starting value effect of ADC. Default is 10.

4. 'SCALE'
<br />Global scaling factor for the input signal. Default is 1. Since here we only demostrate 1 bit quantizer ADC, there is not too much effect.

5. 'BITS'
<br />How many bits of quantizer for ADC. Here we only demonstrate 1 bit quantizer ADC so default is 1. It will have not too much effect untill in the future development multi-bit quantizer is used.

7. 'FILENAME'
<br />The file name in which you want to save or load the random matrix PHI. This file name works for either save and load mode.

8. 'SAVE'
<br />Indicate you want to save the random matrix PHI in 'FILENAME'.

9. 'LOAD'
<br />Indicate you want to load the random matrix PHI from file in 'FILENAME'. This is used to recover previous running results.

10. 'binFFT'
<br />MATLAB vector at which frequency the sparse signals are. Note this is relative to the sampling frequency. For example, default binFFT is [3 51] and default DATA_LENGTH is 256. So actually the two sine signals at frequency of 3\*1/256 and 51\*1/256 are superpositioned as the input signal X(t). The length of binFFT is the sparsity of input signal. It is eqault to the 'S' we mentioend above and in the theoretical background. Note that in addition to satisfaction of inequality <img src="http://latex.codecogs.com/svg.latex?m\geq%20C\times%20S\log%20n" border="0" align="center"/>, different signal's frequency are recommended to be far away from each other. For example, [3 51 97] has much higher probability to be recovered than [3 4 5]. Note that due to [nyquist sampling theorem](https://en.wikipedia.org/wiki/Nyquist%E2%80%93Shannon_sampling_theorem), the maximum elment of 'binFFT' cannot exceeds DATA_LENGHT/2.

11. 'amp'
<br />The amplitude which each signal of binFFT is at. Default is [0.8 0.3] for binFFT of [3 51]. It should has the same length of binFFT.

12. 'alg_choice'
<br />Must be either 1 or 2 or 3. Default is 1 and SPGL1 will be used. 2 is for OMP and 3 is for CoSaMP.

13. 'bp_verbosity'
<br />Verbosity for SPGL1 to ouput log or not. 0=quiet, 1=some output, 2=more output. Default is 1.

14. 'bp_iterations'
<br />Max number of iterations SPGL1 is allowed. Default is 1e07 which is large enough to find a solution in most cases. Usually the SPGL1 will converge in 1000 to 3000 steps.

15. 'bp_nPrevVals'
<br />Number previous func values for linesearch in SPGL1. Default is 3.

16. 'bp_bpTol'
<br />Tolerance for basis pursuit solution in SPGL1. Default is 1e-8 and is accurate enough in most cases.

17. 'bp_optTol'
<br />Optimality tolerance in SPGL1. Default is 1e-6

18. 'bp_decTol'
<br />Required relative change in primal objective. Used for Newton methods in SPGL1. Default is 1e-6.

19. 'bp_stepMin'
<br />Minimum spectral step in SPGL1. Default is 1e-16.

20. 'bp_stepMax'
<br />Maximum spectral step in SPGL1. Default is 1e5.

21. 'bp_subspaceMin'
<br />0=no subspace minimization, 1=subspace minimization. Default is 0.

22. 'OMP_residual'
<br />The residual value on which the OMP algorithm will stop. Default is 1e-6. Since OMP is a pure greedy search algorithm, when number of entries in binFFT is too high, OMP may not work.

23. 'CoSaMP_estimate_sparsity'
<br />The estimated sparsity of input signal for CoSaMP. Since the convergence criteria is inverse proportional to the sparsity, usually a higher value than you really want to recover is recommended. For example, as default binFFT is [3 51], default 'CoSaMP_estimate_sparsity' is 5 (the actual value should be 2 instead).

## Results and Discussion
Next we will discuss the results on default settings. In default settings, PATHS is 15, DATA_LENGTH is 256 and binFFT is [3 51]. Note that due to random nature of PHI and recovery algorithm, you may not repeat the exact same results. However, when algorithm cannot fully recover all the signals, do not just give up. Try to run once more and you may see different results!

After running command ```SD_COMP_DEMO(struct())``` (running demo with default settings), three graphs will be generated.<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/Figure_1.png"></p>

Figure 1 shows the FFT graph of the original signal in dB20 scale. Here, frequency 3 and 51 are shown where other small values are simply quantization noise (we call them noise floor).

<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/Figure_2.png"></p>

The upper subplot of figure 2 shows the FFT graph of the recovered signal (in blue) and the original signal (in red) together. It is shown in dB20 scale. Since the two spikes of blue overlap with two red spikes, we are successful to recover the compressed sensing signals!

The lower graph shows the difference between original and recovered signals. Note that the different is less than -4dB. That means the original and recovered signals are very close indeed!

<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/Figure_3.png"></p>

The upper subplot of figure 3 shows the time domain graph of the original signal (in blue) and the recovered signal (in green). Note that they are very close indeed.

The lower subplot of figure 3 shows the time domain graph of difference between original signal and recovered signal. The amplitude of error is within 0.2 whereas the original signal amplitude is around 1. This is a very good recovery considering only one bit quantizer is used for sigma delta ADC.

In conclusion, we recovered a signal of lenght of 256, which is a superposition of sine wave of frequency 3 and 51, with 15 paths of ADC. By transmitting only this 15 paths' accumulated output, we are able to recover the input signal of sparsity of 2. The compression ration is around 15/256=5.86%! That is amazing!

You may think since m is only equal to 2 in this example and the next question is is it really capable of reconstructing more sparse signals? Well, I tried with PATHS=45, DATA_LENGTH=1024, binFFT=[3 51 109 211 379 483]. 

In this reconstruction problem, I find that SPGL1 is working well with option ```opts=struct('PATHS',60,'DATA_LENGTH',1024,'amp',[0.8,0.6,0.5,0.6,0.6,0.6],'binFFT',[3,51,109,211,379,483],'alg_choice',1,'verbosity',1)```

OMP is working well with option ```opts=struct('PATHS',60,'DATA_LENGTH',1024,'amp',[0.8,0.6,0.5,0.6,0.6,0.6],'binFFT',[3,51,109,211,379,483],'alg_choice',2)```

CoSaMP at heart is a greedy search algorithm but it incorporates ideas from the combinational algorithms. Therefore it is working well with option ```opts=struct('PATHS',60,'DATA_LENGTH',1024,'amp',[0.8,0.6,0.5,0.6,0.6,0.6],'binFFT',[3,51,109,211,379,483],'alg_choice',2,'CoSaMP_estimate_sparsity',15)```

With either of the three algorithm, reconstruction should be sucessful and this is what I get:

<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/Figure_1_1.png"></p>

<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/Figure_2_1.png"></p>

<p align="center"><img src="https://github.com/ssw5075839/COMP/blob/master/pics/Figure_3_1.png"></p>

This concludes our discussion of compressed sensing sigma delta ADC.
