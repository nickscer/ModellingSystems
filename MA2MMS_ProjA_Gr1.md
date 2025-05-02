---
geometry: margin=2cm
---
[comment]: # (This is how you can make comments without affecting the output file or preview, the above section is the YAML btw, it is not rendered either but rather sets custom global markdown formatting for pdf/html output)

# MA2MMS Project A. Modelling Biological Systems 
*Nick Scerbacenco, Keisha M. Patel and Henry Davis*
*University of Reading*  


## Abstract
Add a brief abstract with a description and conclusion here ...


## Table of Contents

1. [Description of the Biological/Ecological System](#description-of-the-biologicalecological-system)  
   1.1. [Introduction](#introduction)
   1.2  [Model selection](#introduction)
2. [Stability of the model](#stability-of-the-model)
3. [Appendix](#appendix)  


# Description of the Biological/Ecological System

Add some description 

## Introduction

## Model selection
 
We are looking at ... red and grey squirrel populations over time using the Lotka-Volterra Competitive model. 

$$ \frac{\mathrm{d}R}{\mathrm{d}t} = r_{R}R(1-\frac{R-G\alpha_{RG}}{K_{R}}) $$
$$ \frac{\mathrm{d}G}{\mathrm{d}t} = r_{G}G(1-\frac{G-R\alpha_{GR}}{K_{G}}) $$

where $R(t)$ and $G(t)$ represent the population of red and grey squirrels at a given time $t$,
$r_{R}$ and $r_{G}$ represent the intrinsic growth rates of red and grey squirrels respectively,
$K_{R}$ and $K_{G}$ represent the carrying capcity of red and grey squirrels,
and finally $\alpha_{RG}$ and $\alpha_{GR}$ represent the competition coefients - the effect of grey squirrels on red squirrels and the effect of red squirrels on grey squirrels respectively.

## Stability of the model

$$ \frac{\mathrm{d}R}{\mathrm{d}t} = 0.61R(1-\frac{R-0.8G}{K_{R}}) $$
$$ \frac{\mathrm{d}G}{\mathrm{d}t} = 0.82G(1-\frac{G-0.09R}{K_{G}}) $$

where $K_{G} = 3×10^{6}$ and $K_{R} = 2.5×10^{6}$

Calculating the Jacobian:

$$ J(R,G) = \begin{bmatrix} \frac{\partial{F}}{\partial{R}} & \frac{\partial{F}}{\partial{G}} \\
\frac{\partial{E}}{\partial{R}} & \frac{\partial{E}}{\partial{G}}
\end{bmatrix} $$

Where $F(R,G) = \frac{\mathrm{d}R}{\mathrm{d}t}$ and $E(R,G) = \frac{\mathrm{d}G}{\mathrm{d}t}$.
Substituting in the partial derivatives we have:

$$ J(R,G) = \begin{bmatrix} 0.61 - \frac{1.22}{K_{R}}R - \frac{0.488}{K_{R}}G & -\frac{0.488}{K_{R}} \\
-\frac{0.0738}{K_{G}}G & 0.82 - \frac{0.18}{K_{G}}G - \frac{0.0738}{K_{G}}R \end{bmatrix} $$

we calculated the nullclines by setting $F(R,G) = 0$ and $E(R,G) = 0$ (.. link to working out in appendix) and got the equilibria points to be.. $(0,0)$, $(K_{R},0)$, $(K_{G},0)$ and $(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928})$.

Now looking at the Jacobian at the equilibria points:

$ J(0,0) = \begin{bmatrix} 0.61 & 0 \\
0 & 0.82 \end{bmatrix} $

.Examining the stability by calculating the eigenvalues. Since $J(0,0)$ is a diagonal matrix, $\lambda_{1} = 0.61$ and $\lambda_{2} = 0.82$. Therefore this point is a unstable source;

$$J(K_{R},0) = \begin{bmatrix} -0.61 & -0.488 \\
0 & 0.82 - \frac{0.073K_{R}}{K_{G}}\end{bmatrix} $$
.Examining the stability by calculating the eigenvalues. Since $J(K_{R},0)$ is an upper triangular matrix, $\lambda_{1} = -0.61$ and $\lambda_{2} = 0.82 - \frac{0.073K_{R}}{K_{G}}$. Therefore this point is a saddle point;


# Appendix

*Add any ideas/manuscripts, links and references below, treat this as a draft for now*
