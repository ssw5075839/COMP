#Compressed-Sensing Sigma Delta ADC with SPG, OMP and CoSaMP
This is a MATLAB simulation program for compressed-sensing sigma-delta analog to digital converter (ADC).

##Introduction to Basis Pursuit Denoising

The idea of compressed-sensing comes from [basis pursuit denoising (BPDN)](https://en.wikipedia.org/wiki/Basis_pursuit_denoising) problem in machine learning. Basically, this problem has three equivalent forms:

1. Basis pursuit denoising:
<p align="center">
Minimize <img src="http://www.sciweavers.org/upload/Tex2Img_1484544777/render.png" align="center" border="0" alt="\left \| x \right \|_1" width="39" height="19" /> subject to <img src="http://www.sciweavers.org/upload/Tex2Img_1484545044/render.png" align="center" border="0" alt="\left \| A \overrightarrow{x} - \overrightarrow{y}  \right \|_2 \leq  \sigma " width="135" height="35" />
</p>
2. Basis pursuit:
<p align="center">
Minimize <img src="http://www.sciweavers.org/upload/Tex2Img_1484545104/render.png" align="center" border="0" alt="\left \| x \right \|_1" width="39" height="19" /> subject to <img src="http://www.sciweavers.org/upload/Tex2Img_1484545131/render.png" align="center" border="0" alt="A \overrightarrow{x} = \overrightarrow{y} " width="75" height="26" />
</p>
3. Lasso:
<p align="center">
Minimize <img src="http://www.sciweavers.org/upload/Tex2Img_1484545155/render.png" align="center" border="0" alt="\left \| A \overrightarrow{x} - \overrightarrow{y}  \right \|_2" width="100" height="35" /> subject to <img src="http://www.sciweavers.org/upload/Tex2Img_1484545178/render.png" align="center" border="0" alt="\left \| x  \right \|_1 \leq  \tau " width="74" height="19" />
</p>

Machine learning researchers have developed various kinds of algorithms to solve BPDN problem. Some basic textbook materials can be found in Kevin P. Murphy's [Machine Learning: A Probabilistic Perspective](https://www.amazon.com/Machine-Learning-Probabilistic-Perspective-Computation/dp/0262018020/ref=sr_1_1?ie=UTF8&qid=1484504311&sr=8-1&keywords=machine+learning+a+probabilistic+perspective) Chapter 13.2.3 sparse linear models.

In this project, we are going to use three algorithms implmente in MATLAB:

1. SPGL1: [Spectrial Projected Gradient for L1 minimization](https://github.com/mpf/spgl1)<br /><br />This algorithm is developed by Professor Michael P. Friedlander fron University of British Columbia. The details can be found in his paper [Probing the Pareto frontier for basis pursuit solutions](https://www.cs.ubc.ca/~mpf/pubs/probing-the-pareto-frontier-for-basis-pursuit-solutions/). As a general comment on this implementation, it is an iterative algorithm that can solve large scale BPDN reconstruction problem. However, even for small scale problem it still runs for many iterations and takes a relative long time to find the root.
2. OMP: Orthogonal Matching Pursuit<br /><br /> This algorithm is proposed by [Stephane G. Mallat] (https://www.di.ens.fr/~mallat/papiers/MallatPursuit93.pdf) in 1993. It is a greedy search algorithm on BPDN least square equation. It is suitable for small scale problem when dimension of X is small and computing requirement is affordable. In this scenario, it is much faster than SPGL1 but we do observe that it has some accuracy problem. Due to its greedy nature, OMP sometimes could not find all the features of X while SPGL1 could. In this project, we use [Stephen Becker](https://www.mathworks.com/matlabcentral/fileexchange/32402-cosamp-and-omp-for-sparse-recovery)'s implementation in MATLAB.
3. CoSaMP:Compressive Sampling Matched Pursuit<br /><br />This algorithm is from [D. Needell, J. A. Tropp](https://arxiv.org/pdf/0803.2392v2.pdf). At hearts, this algorithm is greedy pursuits but it also incorporates ideas from the combinatorial algorithms to guarantee speed and to provide rigorous error bounds. The details can be found in the paper. In this project, we use [Stephen Becker](https://www.mathworks.com/matlabcentral/fileexchange/32402-cosamp-and-omp-for-sparse-recovery)'s implementation in MATLAB.

##Introduction to Compressed Sampling
Here, I will present a highly simplified introduction to the compressed sampling. For more details and rigorous mathematical proof, please refer to the paper [An Introduction To Compressive Sampling](http://dsp.rice.edu/sites/dsp.rice.edu/files/cs/CSintro.pdf).

1. What is signal sensing or sampling:<br /><br />In this artical, we will discuss sensing mechanisms in which information about a time domain signal f(t) is obtained by a series of recording values <p align="center"> <img src="http://www.sciweavers.org/upload/Tex2Img_1484592384/render.png" align="center" border="0" width="114" height="19" /> </p> In another word, we correlate signal f(t) with a series of pre-defined standard signals <img src="http://www.sciweavers.org/upload/Tex2Img_1484582992/render.png" align="center" border="0" alt="\varphi_{k}(t)" width="44" height="19" />, for which we have various kinds of choices. For example, if <img src="http://www.sciweavers.org/upload/Tex2Img_1484582992/render.png" align="center" border="0" alt="\varphi_{k}(t)" width="44" height="19" /> are Dirac delta functions (spikes), then we are simpling recording f(t) with its discrete sampled values y. If <img src="http://www.sciweavers.org/upload/Tex2Img_1484582992/render.png" align="center" border="0" alt="\varphi_{k}(t)" width="44" height="19" /> are sine wave functions <img src="http://www.sciweavers.org/upload/Tex2Img_1484583389/render.png" align="center" border="0" alt="sin( \omega _k t)" width="72" height="18" /> in which <img src="http://www.sciweavers.org/upload/Tex2Img_1484583501/render.png" align="center" border="0" alt="\omega_k = k \times  \omega " width="93" height="18" />, we are actually transforming f(t) into frequency domain and this process is called discrete Fourier transformation. We give these pre-defined standard signals <img src="http://www.sciweavers.org/upload/Tex2Img_1484582992/render.png" align="center" border="0" alt="\varphi_{k}(t)" width="44" height="19" /> a new name: basis function. If they are orthonormal basis, we call them orthobasis. The sine wave functions described above is an orthobasis example.
2. Sparse signal and sparsity:<br /><br />Mathematically speaking, in previous section, we have a signal f in t-domain (note this domain is not limited to time domain and can be extended to other domain, for example, two-dimension space domain for images) and we expand f in an orthonormal basis <img src="http://www.sciweavers.org/upload/Tex2Img_1484585153/render.png" align="center" border="0" alt=" \Phi  = [ \varphi _1,  \varphi _2, ... ,   \varphi _3 ]" width="156" height="19" /> as follows:<p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484592153/render.png" align="center" border="0" alt="f(t) =  \sum_{i=1} ^ n x_i   \varphi  _i(t)" width="133" height="50" /></p>By this expansion, we transform f from t-domain to a new domain and represent it with a series of coefficient <img src="http://www.sciweavers.org/upload/Tex2Img_1484585499/render.png"align="center" border="0" width="106" height="19" />. f(t) may have vary complex form in t-domain, however, it could be sparse in the new domain. That is, many entry of <img src="http://www.sciweavers.org/upload/Tex2Img_1484585806/render.png" align="center" border="0" alt="x_i" width="18" height="15" /> is actually 0 or very close to 0 and only S entries are nonzero. We denote these nonzero entries in a set X<sub>s</sub> and call S as the sparsity of the original signal f(t). To make it more clear, we consider <img src="http://www.sciweavers.org/upload/Tex2Img_1484591388/render.png" align="center" border="0" alt="f_s(t)" width="39" height="19" /> obtained by keeping only S nonzero terms. Here, we define<p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484591858/render.png" align="center" border="0" alt="f_s: = \Psi X_S" width="79" height="19" /></p>where from here, X<sub>s</sub> is the vector of  coefficients (X<sub>i</sub>) with all but S entries set to zero. We say f<sub>S</sub> is the sparse and approximate representation of the original signal f(t) and we have <p align="center"> <img src="http://www.sciweavers.org/upload/Tex2Img_1484592716/render.png" align="center" border="0" alt="\parallel f-f_s \parallel _{l_2} = \parallel X-X_s \parallel _{l_2}" width="182" height="19" /></p>If coefficient X is sparse and compressible in the sense that the sorted magnitude of X<sub>i</sub> decay quickly, then X should be well approximated by X<sub>x</sub> and the error <img src="http://www.sciweavers.org/upload/Tex2Img_1484593014/render.png" align="center" border="0" alt="\parallel f-f_s \parallel _{l_2}" width="81" height="19" /> should be small.<br /><br />If you are farmilar with image processing, JPEG format is a very good example of this principle.
3. Orthobases coherence:<br /><br />From previous section, we know that for a given signal f, it can have different representations in different orthobases. Let us consider two orthobases <img src="http://www.sciweavers.org/upload/Tex2Img_1484594785/render.png" align="center" border="0" alt="( \Phi , \Psi )" width="50" height="18" />, we define the coherence between these two orthobases as <p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484595114/render.png" align="center" border="0" width="272" height="33" /></p>Basically, the coherence is checking how correlated those two orthobases are. You can image that if two orthobases are the same, the coherence is the maximum possible <img src="http://www.sciweavers.org/upload/Tex2Img_1484595351/render.png" align="center" border="0" alt="\sqrt{n}" width="31" height="19" />.<br /><br />How about the least coherent orthobases? One example is the time domain spike basis <img src="http://www.sciweavers.org/upload/Tex2Img_1484595494/render.png" align="center" border="0" alt=" \varphi _k= \delta (t-k)" width="106" height="19" />. This is the most common basis we use in real life. When you record a signal's waveform, you are actually using the spike basis. Its reciprocal basis is Fourier basis, <img src="http://www.sciweavers.org/upload/Tex2Img_1484595669/render.png" align="center" border="0" alt=" \psi _j(t) =  \sqrt{n} e^{i 2 \pi  jt/n}" width="144" height="25" />, has coherence <img src="http://www.sciweavers.org/upload/Tex2Img_1484595750/render.png" align="center" border="0" alt=" \mu ( \Phi , \Psi )=1" width="94" height="19" />. This is the smallest coherence we can have and thus we say spike basis and Fourier basis are maximally incoherent.<br /><br />Finally, we would like point out that a random matrices are largely incoherent with any fixed basis <img src="http://www.sciweavers.org/upload/Tex2Img_1484595998/render.png" align="center" border="0" alt=" \Psi " width="15" height="14" />. We will omit the rigorous proof here (if you are interested, you can refer to the paper I mentioned at the beginning) but just stated: with high probability, an orthobasis <img src="http://www.sciweavers.org/upload/Tex2Img_1484596173/render.png" align="center" border="0" alt="\Phi " width="14" height="14" /> uniformly at random is highly incoherent with any fixed basis <img src="http://www.sciweavers.org/upload/Tex2Img_1484595998/render.png" align="center" border="0" alt=" \Psi " width="15" height="14" />. And the coherence is about <img src="http://www.sciweavers.org/upload/Tex2Img_1484596346/render.png" align="center" border="0" alt=" \mu ( \Phi , \Psi ) \approx  \sqrt{2logn} " width="149" height="26" />. As an exmaple, we will use a matrix <p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484596664/render.png" align="center" border="0" alt="A = \begin{bmatrix}a_{11} & ... & a_{1n} \\... & ...\\a_{m1} & ... & a_{mn} \end{bmatrix} " width="172" height="65" /></p> with entry a<sub>ij</sub> equal to +1 or -1 randomly. This matrices is highly incoherent with spike basis or Fourier basis.
4. Sparse signal reconstruction:<br /><br />Ideally, we would like to collect or measure all the n coefficient of f in our observation orthobasis. However, if we only collect part of the data<p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484606339/render.png" align="center" border="0" width="186" height="43" /></p>M is a subset of cardinality <img src="http://www.sciweavers.org/upload/Tex2Img_1484606474/render.png" align="center" border="0" alt="m < n" width="54" height="12" />.Here, the only other information we know is in some domain <img src="http://www.sciweavers.org/upload/Tex2Img_1484606564/render.png" align="center" border="0" alt=" \Psi " width="15" height="14" />the signal f has a sparse representation x. Maybe with the partial observation, we cannot reconstruct the exact signal f in general. But we can find the approximate reconstrution <img src="http://www.sciweavers.org/upload/Tex2Img_1484606832/render.png" align="center" border="0" alt="f^*" width="22" height="21" /> given by <img src="http://www.sciweavers.org/upload/Tex2Img_1484606871/render.png" align="center" border="0" alt="f^*= \Psi x^*" width="78" height="21" />, where <img src="http://www.sciweavers.org/upload/Tex2Img_1484606917/render.png" align="center" border="0" alt="x^*" width="25" height="17" /> is the solution to the optimization problem:<p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484607295/render.png" align="center" border="0" width="431" height="21" /></p>By minimize the L0 norm of <img src="http://www.sciweavers.org/upload/Tex2Img_1484607403/render.png" align="center" border="0" alt=" \widetilde{x} " width="15" height="17" />, we actually mean to find as sparse solution as possible. However, we know that minimization of L0 norm is not a convex problem, instead, minimization of L1 norm is. Furthermore, we do know that minimization of L1 norm do propose a sparse solution and it is equivalent to mminimization of L0 norm. Therefore, with a single change from L0 to L1, we have a solvable convex optimization problem:<p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484607684/render.png" align="center" border="0" width="425" height="21" /></p> Next, we will prove that if f is sufficiently sparse, the recovery via L1 norm minimization is actually exact rather than apporximate.
5. Theorem 1:<br /><br />Suppose the coefficient vector x of f in the basis <img src="http://www.sciweavers.org/upload/Tex2Img_1484595998/render.png" align="center" border="0" alt=" \Psi " width="15" height="14" /> is S-sparse. Select m measurement in the observation domain <img src="http://www.sciweavers.org/upload/Tex2Img_1484608533/render.png" align="center" border="0" alt=" \phi " width="17" height="19" /> and <img src="http://www.sciweavers.org/upload/Tex2Img_1484608533/render.png" align="center" border="0" alt=" \phi " width="17" height="19" /> is uniformly random. Then if <p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484608687/render.png" align="center" border="0" alt="m \geq C \times  \mu ^2( \Phi , \Psi ) \times S \times \log n" width="232" height="22" /></p> for some positive constant C, the solution in the previous section gives the exact reconstrcution of f with overwhelming probability. The probability will exceeds <img src="http://www.sciweavers.org/upload/Tex2Img_1484608914/render.png" align="center" border="0" alt="1- \delta " width="44" height="17" /> if <img src="http://www.sciweavers.org/upload/Tex2Img_1484608866/render.png" align="center" border="0" alt="m \geq C \times  \mu ^2( \Phi , \Psi ) \times S \times \log (n/ \delta)" width="258" height="22" />. We will skip the proof here. For details, please refer to the paper I mentioned. Here, I interpret the inequality in three points:<br /><br />First, the coherence between two bases are very important. To find a relative incoherent basis with respect to any fixed basis, a random matrix <img src="http://www.sciweavers.org/upload/Tex2Img_1484608533/render.png" align="center" border="0" alt=" \phi " width="17" height="19" /> should be selected. However, due to the fact that machine generated random number is not really random, the quality of randomness will determine the reconstruction quality later on.<br /><br />Second, when the coherence <img src="http://www.sciweavers.org/upload/Tex2Img_1484595750/render.png" align="center" border="0" alt=" \mu ( \Phi , \Psi )=1" width="94" height="19" /> is close to 1, then on the order of <img src="http://www.sciweavers.org/upload/Tex2Img_1484609519/render.png" align="center" border="0" alt="S \log n" width="58" height="19" /> smaples suffice to reconstrut instead of n in normal sensing. The data storage space or transmitting bandwidth save is <img src="http://www.sciweavers.org/upload/Tex2Img_1484609880/render.png" align="center" border="0" alt="S  \frac{1}{n} \log n " width="72" height="43" />.<br /><br />Finally, the only priori we know is f is sparse in domain <img src="http://www.sciweavers.org/upload/Tex2Img_1484595998/render.png" align="center" border="0" alt=" \Psi " width="15" height="14" />. We don't assume any other knowledge about the exact number of S, or their locations, or their amplitude. We simply run L1 minimization algorithm to find a sparse solution. And if the signal is truely sparse enough, the exact reconstruction occurs.
6. Denoising problem:<br /><br />In the previous five sections, we examined the reconstruction without noise. However, in real world, every kind of measurement or sensing introduces noise. Hence, will consider the previous problem in a more generous form.<br /><br />Consider <img src="http://www.sciweavers.org/upload/Tex2Img_1484610453/render.png" align="center" border="0" alt="f =  \Psi x" width="57" height="19" /> and f is S-sparse in domain <img src="http://www.sciweavers.org/upload/Tex2Img_1484595998/render.png" align="center" border="0" alt=" \Psi " width="15" height="14" />. We get the random sensing of the signal f in form of <img src="http://www.sciweavers.org/upload/Tex2Img_1484610593/render.png" align="center" border="0" alt="y = R \Phi f" width="68" height="19" /> where R is a m-by-n random matrix that extract data randomly. Then without noise, we can write <p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484610718/render.png" align="center" border="0" alt="y = Ax, \ where A=R \Phi  \Psi " width="197" height="19" /></p>With noise, we can write <p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484610823/render.png" align="center" border="0" alt="y = Ax+z, \ where A=R \Phi  \Psi \ and \ z \ is \ noise" width="351" height="19" /></p>If R is identity matrix, namely, we take the exact measurment, then A is isometric and is not hard to prove that <img src="http://www.sciweavers.org/upload/Tex2Img_1484611122/render.png" align="center" border="0" alt=" \parallel x \parallel _{l_2} =  \parallel Ax \parallel  _{l_2}" width="124" height="19" />. However, if R is not identity, A is not rigorous isometric and we need to define restricted isometry property (RIP) of A.<br /><br />For each integer S=1,2,..., define the isometry constant <img src="http://www.sciweavers.org/upload/Tex2Img_1484611308/render.png" align="center" border="0" alt=" \delta _s" width="21" height="18" /> of matrix A as the smallest number such that:<p align="center"><img src="http://www.sciweavers.org/upload/Tex2Img_1484611428/render.png" align="center" border="0" alt="(1- \delta _s) \parallel x \parallel _{l_2}^2 \leq  \parallel Ax \parallel _{l_2}^2 \leq (1+ \delta _s) \parallel x \parallel _{l_2}^2" width="321" height="28" /> </p>holds for all the S-sparse vectors x. We will loosely say that a matrix A obeys RIP of order S if <img src="http://www.sciweavers.org/upload/Tex2Img_1484611308/render.png" align="center" border="0" alt=" \delta _s" width="21" height="18" /> is not too close to one. With RIP, A can approximately preserves the Euclidean length of S-sparse vector x. In turn, this implies that x cannot be in the null space of A. Otherwise, our reconstruction of x is impossible.<br /><br />Then, how does RIP relate to our reconstruction problem? Imagine that we get a measurment result y which in truth is equal to Ax<sub>1</sub>. Then we find another root x<sub>2</sub> also satisfy y = Ax<sub>2</sub>. Suppose <img src="http://www.sciweavers.org/upload/Tex2Img_1484612393/render.png" align="center" border="0" alt="\delta _{2s}" width="29" height="18" /> is also sufficiently less than one. Then we will have:<p align="center"> <img src="http://www.sciweavers.org/upload/Tex2Img_1484612521/render.png" align="center" border="0" alt="(1- \delta _{2s}) \parallel x_1-x_2 \parallel _{l_2}^2 \leq  \parallel Ax_1-Ax_2 \parallel _{l_2}^2 \leq (1+ \delta _{2s}) \parallel x_1-x_2 \parallel _{l_2}^2" width="494" height="28" /> </p>This means, we can guarantee x<sub>2</sub> is very close to the underlying truth x<sub>1</sub>. Even with noise presented, their different is only on the order of noise. We will formally write the theorem as below:<br /><br />Theorem 2:<br />If RIP of A holds, then the solution <img src="http://www.sciweavers.org/upload/Tex2Img_1484613127/render.png" align="center" border="0" alt="x_*" width="25" height="14" /> to the follwing problem (basis pursuit):<p align="center">Minimize <img src="http://www.sciweavers.org/upload/Tex2Img_1484612989/render.png" align="center" border="0" alt="\left \|  \widetilde{x}  \right \|_1" width="39" height="25" /> subject to <img src="http://www.sciweavers.org/upload/Tex2Img_1484613065/render.png" align="center" border="0" alt="A \widetilde{x} = \overrightarrow{y} " width="67" height="26" /></p> gives a very close approximation of the underlying truth x. If <img src="http://www.sciweavers.org/upload/Tex2Img_1484613262/render.png" align="center" border="0" alt=" \delta _{2s} <  \sqrt{2} -1" width="108" height="26" />, then the solution <img src="http://www.sciweavers.org/upload/Tex2Img_1484613127/render.png" align="center" border="0" alt="x_*" width="25" height="14" /> obeys <p align="center"> <img src="http://www.sciweavers.org/upload/Tex2Img_1484613446/render.png" align="center" border="0" alt=" \parallel x^*-x \parallel _{l_2} \leq C_0 \parallel x - x_S \parallel _{l_1}/ \sqrt{s} \ and  \parallel x^*-x \parallel _{l_1} \leq C_0 \parallel x - x_S \parallel _{l_1}" width="519" height="21" /> </p> for some constant C<sub>0</sub>. When x is actually S-sparse, x is equal to x<sub>s</sub> and the solution <img src="http://www.sciweavers.org/upload/Tex2Img_1484613127/render.png" align="center" border="0" alt="x_*" width="25" height="14" /> is equal to x. Hence, we get the exact construction. Please note that, compare to theorem 1, theorem 2 involves no probability. It is a deterministic theorem.<br /><br />With noise presented, we introduce the third theorem:<br /><br />Theorem 3:<br />If RIP of A holds, then the solution <img src="http://www.sciweavers.org/upload/Tex2Img_1484613127/render.png" align="center" border="0" alt="x_*" width="25" height="14" /> to the follwing problem (basis pursuit denoising):<p align="center">Minimize <img src="http://www.sciweavers.org/upload/Tex2Img_1484612989/render.png" align="center" border="0" alt="\left \|  \widetilde{x}  \right \|_1" width="39" height="25" /> subject to <img src="http://www.sciweavers.org/upload/Tex2Img_1484614197/render.png" align="center" border="0" alt=" \parallel A \widetilde{x} - \overrightarrow{y} \parallel _{l_2} \leq  \epsilon " width="132" height="26" /></p> gives a very close approximation of the underlying truth x. Here, <img src="http://www.sciweavers.org/upload/Tex2Img_1484614287/render.png" align="center" border="0" alt=" \epsilon " width="15" height="12" /> bounds the noise in the sensing of data. If <img src="http://www.sciweavers.org/upload/Tex2Img_1484613262/render.png" align="center" border="0" alt=" \delta _{2s} <  \sqrt{2} -1" width="108" height="26" />, then the solution <img src="http://www.sciweavers.org/upload/Tex2Img_1484613127/render.png" align="center" border="0" alt="x_*" width="25" height="14" /> obeys <p align="center"> <img src="http://www.sciweavers.org/upload/Tex2Img_1484614044/render.png" align="center" border="0" alt=" \parallel x^*-x \parallel _{l_2} \leq C_0 \parallel x - x_S \parallel _{l_1}/ \sqrt{s}+C_1 \times  \epsilon" width="332" height="21" /> </p> for some constant C<sub>0</sub> and C<sub>1</sub>. When x is actually S-sparse, x is equal to x<sub>s</sub> and the solution <img src="http://www.sciweavers.org/upload/Tex2Img_1484613127/render.png" align="center" border="0" alt="x_*" width="25" height="14" /> is equal to x. Hence, we get the exact reconstruction.<br /><br />All the theorms involves long mathematical proof and we omit them for the sake of readability. If you are really interested in the provment of those inequality, please refer to the paper I mentioned above.
