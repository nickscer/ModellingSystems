---
geometry: margin=2cm
---
[comment]: # (This is how you can make comments without affecting the output file or preview, the above section is the YAML btw, it is not rendered either but rather sets custom global markdown formatting for pdf/html output)
[comment]: # (Use `pandoc test.txt -o test.pdf` in the zsh terminal for compilation)

# MA2MMS Project A. Modelling Ecological Systems 

*Nick Scerbacenco, Keisha Markey Patel and Henry Davis*

*University of Reading*  

## Abstract

Add a brief abstract with a description and conclusion here ...

## Table of Contents

1. [Description of the Ecological System](#description-of-the-ecological-system)  
   2. [Introduction](#introduction)
   3. [Model selection](#model-selection)
	   1. [Choosing Coefficients](#choosing-coefficients)
4. [Stability of the model](#stability-of-the-model)
5. [Numerical Solution](#numerical-solution)
	1. [Local Truncation Error](#local-truncation-error)
	2. [Potential Bifurcation](#potential-bifurcation)
6. [Appendix](#appendix)  
	1. [Bibliography](#bibliography) 

# Description of the Ecological System

This report outlines the relationship between the red squirrel and the grey squirrel in the UK. The two species have a very tense relationship due to competition and displacement which has affected their individual populations. We will look at the different numbers of each and format solvable equations for the each of the populations of the two species. As time has gone on, the population of grey squirrels has increased while the population of red squirrels has decreased. First, we need to look at the kind of relationship that the pair have and the varied biotic and abiotic factors that affect their population sizes. While red and grey squirrels do not directly compete against each other, they compete in terms of food resources and habitats. 

## Introduction

Grey squirrels were introduced to the UK in 19th century and over time have outcompeted red squirrels in their habitats. This is due to their larger builds, their adaptability, and immunity from certain diseases, meaning that the grey squirrel can outcompete the red squirrel without direct harm. The grey squirrel primarily outcompetes the red squirrel through the monopolisation of food resources [^2] (Wildlife Online, n.d.), leaving the red squirrel with a limited supply. It is undoubtable that the grey squirrel has populated much of the urban environment while the red squirrel tends to reside in more rural areas. This is a good representation of the grey squirrels’ adaptability, their ability to fit into the environment that they find and even go to the extent of consuming the different foods available there [^2] (Wildlife Online, n.d.), including that of a more urban area. 

Grey squirrels currently populate the UK in many places, carrying diseases such as squirrel pox, which although is not very harmful to them, it is fatal if passed on to a red squirrel [^1] (Red Squirrel Survival Trust, n.d.). This, in turn, decreases the number of red squirrels while the population of the grey squirrels is only marginally decreased. Alongside this, they are significantly larger than the red squirrel, making them the more likely survivor when in combat with competitors or predators, such as the pine marten. Due to the red squirrel's sensitive immune system, they have on average a lower life expectancy than that of the red squirrel. The red squirrel typically lives up to three years [^1] (Red Squirrel Survival Trust, n.d.), while the grey squirrel usually lives up to nine years in the wild [^4] (Wytham Woods, n.d.). Both species have been known to live up to ten or even fifteen years old, however, the contributions of the outside world, especially the different diseases they face, lowers this significantly. 

Currently, the grey squirrel is the dominant species, while the red squirrel population is declining. However, the population of the red squirrel is slowly rising as the population of the pine marten increases. The pine marten is the natural predator to the grey squirrel, which has a “strong negative” [^3] (British Red Squirrel, n.d.) impact on their numbers, allowing the red squirrel to benefit and thrive as there is less competition for the red squirrel. After being pushed to near extinction, the red squirrel finds peaceful habitat in places less abundant with grey squirrels, meaning their numbers can slowly increase. Evidence was found to suggest that rather than the grey squirrel being hunted by the pine marten to decrease their numbers, the grey squirrel started to migrate to areas less populated with the pine marten. The red squirrel, in turn, benefitted from this, as they were given more land and area to populate without being disturbed by the grey squirrel and their numbers stared to increase. Furthermore, the culling of grey squirrels in legal in the UK as a way of controlling their population size and allowing the population of the red squirrel to increase [^5] (British Red Squirrel, n.d.). This includes trapping, shooting, fertility control and more. 

Looking into the history of the introduction of both of the species, we can start to compare the numbers and substitute them into the Lotka-Volterra Competition model. The introduction of the grey squirrel in England was in 1876 and continued to be released until the 1920’s when the destruction of the red squirrel population was noticed and it soon became illegal to release a grey squirrel into the wild [^6] (British Red Squirrel, n.d.). Red squirrels, however, are native to the UK and lived in the UK prior to the grey squirrel for around 10,000 years [^7] (The Wildlife Trusts, 2020). It is estimated that at the time of the introduction of grey squirrels, the red squirrel population stood at roughly 3,500,000 and has decreased to 140,000 in the past years, while the grey squirrel population has increased to approximately 2,520,000 recorded in 2009 [^8] (Aebischer, Davey and Kingdon, 2011). We have taken the carrying capacity of grey squirrels to be 2,500,000 and used this to find that the carrying capacity for red squirrels is 3,000,000 [^7] (The Wildlife Trusts, 2020).



## Model Selection 

Our aim is modelling *global* populations of red and grey squirrels over time. We begin by considering the *generic* *Lotka-Volterra* system s.t. 

$$
\overset{\text{Generic Model}}{
\boxed{
\begin{aligned}
\frac{\mathrm{d}x}{\mathrm{d}t} &= x(r_1 + a_{11}x + a_{12}y), \\
\frac{\mathrm{d}y}{\mathrm{d}t} &= y(r_2 + a_{21}x + a_{22}y).
\end{aligned}
}
}
$$

where $x$ and $y$ would represent the populations of red and grey squirrels respectively,

- $r_{1}$ - idealistic growth rate of red squirrels 
- $r_{2}$ - idealistic growth rate of grey squirrels 
- $a_{11}$ - measure of limitation on red squirrels
- $a_{22}$ - measure of limitation on grey squirrels 
- $a_{12}$ - measure of competition towards red squirrels i.e. how much $y$s hurt $x$s 
- $a_{21}$ - measure of competition towards greys squirrels i.e. how much $x$s hurt $y$s

${N}\mkern -8.2mu\textcolor{red}{{B}}$ $a_{11}$ and $a_{22}$ along with $r_{1}$ and $r_{2}$ respectively yield the carrying capacities $K_{R}$ and $K_{G}$.

The model in this form is very clear and versatile as it allows easy implementation into numerical algorithms because the equations do not involve division.
 
${N}\mkern -8.2mu\textcolor{red}{{B}}$ This form will prove to be convenient when we will solve the system numerically later on.

For now we introduce a new form of our model with the carrying capacities $K_{R}$ and $K_{G}$ explicitly present to aid analysis of stability and picking coefficients below.      

### Choosing Coefficients
 
We are looking at ... red and grey squirrel populations over time using the *generic* *Lotka-Volterra* Competitive model in a *refined form*. 

$$
\overset{\text{Generic Model - Refined form}}{
\boxed{
\begin{aligned}
\frac{\mathrm{d}R}{\mathrm{d}t} &= r_{R}R\left(1-\frac{R+\alpha_{RG}G}{K_{R}}\right) \\
\frac{\mathrm{d}G}{\mathrm{d}t} &= r_{G}G\left(1-\frac{G+\alpha_{GR}R}{K_{G}}\right) 
\end{aligned}
}}
$$

where $R(t)$ and $G(t)$ represent the population of red and grey squirrels at a given time $t$,
$r_{R}$ and $r_{G}$ represent the intrinsic growth rates of red and grey squirrels respectively,
$K_{R}$ and $K_{G}$ represent the carrying capcity of red and grey squirrels,
and finally $\alpha_{RG}$ and $\alpha_{GR}$ represent the competition coefients - the effect of grey squirrels on red squirrels and the effect of red squirrels on grey squirrels respectively.

${N}\mkern -8.2mu\textcolor{red}{{B}}$ This model is based on the logistic growth model ($\frac{\mathrm{d}x}{\mathrm{d}t} = rx ( 1 - \frac{x}{K} )$ ), with the addition of competition between the two species as they compete for the same natural resources. 

*add more reasoning and sources behind coeffs*

## Stability of the model

We have the coupled differential equations

$$ 
\frac{\mathrm{d}R}{\mathrm{d}t} = 0.61R(1-\frac{R+0.8G}{K_{R}}) 
$$
$$ 
\frac{\mathrm{d}G}{\mathrm{d}t} = 0.82G(1-\frac{G+0.09R}{K_{G}}) 
$$

where $K_{G} = 3×10^{6}$ and $K_{R} = 2.5×10^{6}$

To find the stability of the model, we have to find the equilibria of the system of equations and examine the stability of these points:

The R-nullclines are found to be $R = 0$ or $R = K_{R} - 0.8G$ by setting $\frac{\mathrm{d}R}{\mathrm{d}t} = 0$.

The G-nullclines are found to be $G = 0$ or $G = K_{G} - 0.09R$ by setting $\frac{\mathrm{d}G}{\mathrm{d}t} = 0$.

The equilibrium points are found from the intersections of the R and G nullclines, and are given to be: $(0,0)$, $(K_{R},0)$, $(K_{G},0)$ and $(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928})$.

The Jacobian matrix of the system is given by:

$$ 
J(R,G) = \begin{bmatrix} 0.61 - \frac{1.22}{K_{R}}R - \frac{0.488}{K_{R}}G & -\frac{0.488}{K_{R}} \\
-\frac{0.0738}{K_{G}}G & 0.82 - \frac{0.18}{K_{G}}G - \frac{0.0738}{K_{G}}R \end{bmatrix} 
$$

Now looking at the Jacobian at the equilibria points:

At the equilibrium point $(0,0)$, we have the Jacobian matrix

$$ 
J(0,0) = \begin{bmatrix} 0.61 & 0 \\
0 & 0.82 \end{bmatrix} 
$$

Since $J(0,0)$ is a diagonal matrix, it has two eigenvalues which are $\lambda_{1} = 0.61$ and $\lambda_{2} = 0.82$. Since these eigenvalues are non-positive real numbers, the corresponding fixed point is an unstable source.
this means....

At the equilibrium point $(K_{R},0)$, we have the Jacobian matrix

$$ 
J(K_{R},0) = \begin{bmatrix} -0.61 & -0.488 \\
0 & 0.82 - \frac{0.073K_{R}}{K_{G}}\end{bmatrix} 
$$

Since $J(K_{R},0)$ is an upper triangular matrix, it has two eigenvalues which are $\lambda_{1} = -0.61$ and $\lambda_{2} = 0.82 - \frac{0.073K_{R}}{K_{G}}$. Since these eignevalues are real numbers with $\lambda_{1} < 0 < \lambda_{2}$, the corresponding fixed point is an unstable saddlepoint.
This means...

At the equilibrium point $(0,K_{G})$, we have the Jacobian matrix

$$ 
J(0,K_{G}) = \begin{bmatrix} 0.61 - \frac{0.488K_{G}}{K_{R}} & 0 \\
-0.0738 & -0.82 \end{bmatrix} 
$$

Since $J(0,K_{G})$ is a lower triangular matrix, it has two eigenvalues which are $\lambda_{1} = 0.61 - \frac{0.488K_{G}}{K_{R}}$ and $\lambda_{2} = -0.82$. Since these eignevalues are real numbers with $\lambda_{2} < 0 < \lambda_{1}$, the corresponding fixed point is an unstable saddlepoint.
This means...

Finally at the equilibrium point $(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928})$, we have the Jacobian matrix

$$ 
J(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928}) = \begin{bmatrix} \frac{0.488K_{G} - 0.61K_{R}}{0.928K_{R}} & \frac{0.3904K_{G} - 0.488K_{R}}{0.928K_{R}} \\
\frac{0.006642K_{R} - 0.0738K_{G}}{0.928K_{G}} & \frac{0.738K_{R} - 0.82K_{G}}{0.928K_{G}} \end{bmatrix} 
$$

$J(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928})$ has two eigenvalues which are 

$$
\lambda_{1,2} = \frac{1}{2} \left[
\frac{0.488 K_{G}^2 + 0.738 K_{R}^2 - 1.43 K_{R} K_{G}}{0.928 K_{R} K_{G}} \pm \sqrt{ \left( \frac{\text{num1}}{0.928 K_{R} K_{G}} \right)^2 -4\frac{\text{num2}}{(0.928)^2 K_{R} K_{G}}}\right]
$$

${N}\mkern -8.2mu\textcolor{red}{{B}}$ See *[Notes on Stability](#notes-on-stability)* for precise calculations of the numerators inside the square root.

Since these eigenvalues are real, negative numbers, ($\lambda_{1} \approx -0.0186$ and $\lambda_{2} \approx -0.2287$) the corresponding fixed point is an asymptotically stable sink point.
This means...

# Numerical Solution

For the purposes of constructing our numerical algorithm, we return to our generic model in its original form. So we must calculate the coefficients $r_{1},r_{2},a_{11},a_{12},a_{21},a_{22}$ for our populations $x$ and $y$ of red and grey squirrels respectively. 

$$
\begin{aligned}
\frac{\mathrm{d}R}{\mathrm{d}t} &= r_{R}R\left(1-\frac{R+\alpha_{RG}G}{K_{R}}\right) \\
\frac{\mathrm{d}G}{\mathrm{d}t} &= r_{G}G\left(1-\frac{G+\alpha_{GR}R}{K_{G}}\right) 
\end{aligned} \qquad \longrightarrow \qquad
\begin{aligned}
\frac{\mathrm{d}x}{\mathrm{d}t} &= x(r_1 + a_{11}x + a_{12}y), \\
\frac{\mathrm{d}y}{\mathrm{d}t} &= y(r_2 + a_{21}x + a_{22}y).
\end{aligned}
$$

We expand the equations and compare coefficients for $x = R$ and $y=R$ s.t. 

$$
\begin{aligned}
\frac{\mathrm{d}R}{\mathrm{d}t} &= r_R R \left(1 - \frac{R}{K_R} - \frac{G \alpha_{RG}}{K_R}\right) = \underbrace{r_R}_{r_1} R + \underbrace{\left(-\frac{r_R}{K_R}\right)}_{a_{11}} R^2 + \underbrace{\left(-\frac{r_R \alpha_{RG}}{K_R}\right)}_{a_{12}} RG \\
\frac{\mathrm{d}G}{\mathrm{d}t} &= r_G G \left(1 - \frac{G}{K_G} - \frac{R \alpha_{GR}}{K_G}\right) = \underbrace{r_G}_{r_2} G + \underbrace{\left(-\frac{r_G \alpha_{GR}}{K_G}\right)}_{a_{21}} RG + \underbrace{\left(-\frac{r_G}{K_G}\right)}_{a_{22}} G^2
\end{aligned}
$$

${N}\mkern -8.2mu\textcolor{red}{{B}}$ We use coefficients found above, $K_{G}=3\cdot10^6$, $K_{R}=2.5\cdot10^6$, $\alpha_{RG}=0.8$ and $\alpha_{GR}=0.09$. 

$$
\text{By comparison we have}  \qquad \boxed{
\begin{aligned}
r_1 &= r_R \\
r_2 &= r_G \\
a_{11} &= -\frac{r_R}{K_R} \\
a_{12} &= -\frac{r_R \alpha_{RG}}{K_R} \\
a_{21} &= -\frac{r_G \alpha_{GR}}{K_G} \\
a_{22} &= -\frac{r_G}{K_G}
\end{aligned}
} \quad \longrightarrow \quad \boxed{
\begin{aligned}
r_{1} &= 0.61 \\
r_{2} &= 0.82 \\
a_{11} &= -2.44\cdot10^{-7} \\
a_{12} &= -1.952\cdot10^{-7} \\
a_{21} &= -2.46 \cdot 10^{-8} \\
a_{22} &= -2.7\dot{3} \cdot 10^{-7} 
\end{aligned}
}
$$

$\textcolor{red}{!}$ The $a_{22}$ is a recurring decimal marked $\dot{3}$. In the actual implementation, 8 decimal points will be used to minimise the *roundoff error*.  

## Local Truncation Error

add discussion on error inherited by the numerical method

## Potential Bifurcation

add some on bifurcation [^99]


# Appendix

*Add any ideas/manuscripts, links and references below, treat this as a draft for now*

## Notes on Stability 

$$
\begin{aligned}
\text{num1} &= 0.488 K_{G}^2 + 0.738 K_{R}^2 - 1.43 K_{R} K_{G} \\
\text{num2} &= (0.488 K_{G} - 0.61 K_{R})(0.738 K_{R}  - 0.82 K_{G}) \\
& -(0.3904 K_{G} - 0.488 K_{R})(0.006642 K_{R} - 0.0738 K_{G})
\end{aligned}
$$

Jacobian formula

$$
J(R,G) =\begin{bmatrix}
\frac{\partial{F}}{\partial{R}} & \frac{\partial{F}}{\partial{G}} \\
\frac{\partial{E}}{\partial{R}} & \frac{\partial{E}}{\partial{G}}
\end{bmatrix}
$$

Where $F(R,G) = \frac{\mathrm{d}R}{\mathrm{d}t}$ and $E(R,G) = \frac{\mathrm{d}G}{\mathrm{d}t}$

## Notes on Runge-Kutta

An $s$-stage **Runge–Kutta method** approximates $u$ by specifying constants $a_{ij}$, $b_j$ 
and $c_i$ for $i,j = 1,2,\dots,s$ where $s \in \mathbb{N}$.

$$
\overset{\text{Compute the $s$ intermediary steps}}{\boxed{k_i = f\bigl(t_{n-1} + c_i\,h,\;U^{n-1} + h\sum_{j=1}^s a_{ij}\,k_j\bigr) \quad \text{for} \quad i = 1,2,\dots,s}} \qquad \overset{\text{Compute the $n$th approximation}}{\boxed{U^n = U^{n-1} + h\sum_{j=1}^s b_j\,k_j}}
$$

## Python code 

### 4th Order Runge-Kutta method 

```python 
import numpy as np
import matplotlib.pyplot as plt
```
to be added 

## Bibliography

[^99]: Bifurcation theory: https://en.wikipedia.org/wiki/Bifurcation_theory
-	[^1]: Red Squirrel Survival Trust. (n.d.). Red and Grey Squirrels – The differences. Available at: https://www.rsst.org.uk/red-and-grey-squirrels-the-differences/
-	[^2]: WildlifeOnline (n.d.). Squirrels Interaction with Other Species – Decline of the Red | Wildlife Online. Available at: https://www.wildlifeonline.me.uk/animals/article/squirrels-interaction-with-other-species-decline-of-the-red#:~:text=Initially%2C%20it%20was%20thought%20that,i.e.%20they%20monopolise%20food%20resources
-	[^3]: British Red Squirrel (n.d.).  Pine Marten. Available at: http://www.britishredsquirrel.org/grey-squirrels/pine-martin/#:~:text=“Our%20study%20has%20confirmed%20that,martens”%2C%20says%20Dr%20Sheehy
-	[^4]: Wytham Woods (n.d.). Grey squirrel (Sciurus carolinensis). Available at: https://www.wythamwoods.ox.ac.uk/article/grey-squirrel-sciurus-carolinensis#:~:text=Whilst%20Grey%20squirrels%20are%20mostly,lifespan%20is%20approximately%209%20years
-	[^5]: British Red Squirrel (n.d.). Grey Management. Available at: http://www.britishredsquirrel.org/grey-squirrels/grey-control/
-	[^6]: British Red Squirrel (n.d.). Grey Squirrels. Available at: http://www.britishredsquirrel.org/grey-squirrels/#:~:text=Grey%20squirrels%20were%20first%20introduced,grey%20squirrel%20to%20the%20wild
-	[^7]: The Wildlife Trusts (2020). Red Squirrels | The Wildlife Trusts. Available at: https://www.wildlifetrusts.org/red-squirrels#:~:text=Red%20squirrels%20are%20our%20native,a%20wild%20population%20is%201876
-	[^8]: Aebischer, N.J., Davey, P.D. and Kingdon, N.G. (2011) Grey squirrel – Game and Wildlife Conservation Trust. Available at: https://www.gwct.org.uk/research/long-term-monitoring/national-gamebag-census/mammal-bags-comprehensive-overviews/grey-squirrel/

---
#### Draft (to be removed)

$$
\lambda_{1,2} = \frac{1}{2} \left[
\frac{0.488 K_{G}^2 + 0.738 K_{R}^2 - 1.43 K_{R} K_{G}}{0.928 K_{R} K_{G}} \pm \sqrt{ \left( \frac{0.488 K_{G}^2 + 0.738 K_{R}^2 - 1.43 K_{R} K_{G}}{0.928 K_{R} K_{G}} \right)^2 -4\frac{(0.488 K_{G} - 0.61 K_{R})(0.738 K_{R} - 0.82 K_{G})-(0.3904 K_{G} - 0.488 K_{R})(0.006642 K_{R} - 0.0738 K_{G})}{(0.928)^2 K_{R} K_{G}}}\right]
$$

