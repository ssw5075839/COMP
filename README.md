#Compressed-Sensing Sigma Delta ADC with SPG, OMP and CoSaMP
This is a MATLAB simulation program for compressed-sensing sigma-delta analog to digital converter (ADC). For compressed sensing theory, please refer to my file [Theoretical_Background.md](https://github.com/ssw5075839/COMP/blob/master/Theoretical_Background.md). In this Readme file, I will only briefly introduce the concepts and then give the code running instruction.

##Compressed Sensing in Sigma Delta Analog to Digital Converter (ADC)
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

##Code Running Instruction
In this matlab demo of compressed sensing sigma delta ADC, we will implment a one-bit first order sigma delta ADC behaviour model together with compressed sensing. Then we will feed ADC outputs in one of the SPGL1, OMP or CoSaMP to try to reconstruct sparse signal. If the signal is sparse enough in frequency domain and ADC path number is large enough, then all the algorithm should be able to reconstruct the original signal.

To run the code, first define a MATLAB structure opts. (This step is optional. If you don't define opts, the code will use all the default settings.)

To define a struct, use command opts = struct('field1',value1,...,'fieldN',valueN). Each option comes with a ('field', value) pair. The available filed name and suggested value will be discussed in the next section.

Then run the function with command SD_COMP_DEMO(opts) (or simply SD_COMP_DEMO() if you prefer to use default settings). The demo will output the reconstruction algorithm running status in the console and the oupt three graphs to compare the original signal agains reconstructed signal in detail, both in time domain and frequency domain.

##Explanation of Parameters in Program
The available option fields and its suggested values are:
1. 'PATHS'
2. 'DATA_LENGTH'
3. 'WARMUP_CYCLES'
4. 'SCALE'
5. 'BITS'
6. 'FS'
7. 'FILENAME'
8. 'SAVE'
9. 'LOAD'
10. 'binFFT'
11. 'amp'
12. 'alg_choice'
13. 'bp_verbosity'
14. 'bp_iterations'
15. 'bp_nPrevVals'
16. 'bp_bpTol'
17. 'bp_optTol'
18. 'bp_decTol'
19. 'bp_stepMin'
20. 'bp_stepMax'
21. 'bp_subspaceMin'
22. 'OMP_residual'
23. 'CoSaMP_estimate_sparsity'
##Results and Discussion
