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

We have the coupled differential equations

$$ \frac{\mathrm{d}R}{\mathrm{d}t} = 0.61R(1-\frac{R-0.8G}{K_{R}}) $$
$$ \frac{\mathrm{d}G}{\mathrm{d}t} = 0.82G(1-\frac{G-0.09R}{K_{G}}) $$

where $K_{G} = 3×10^{6}$ and $K_{R} = 2.5×10^{6}$

To find the stability of the model, we have to find the equilibria of the system of equations and examine the stability of these points:

The R-nullclines are found to be $R = 0$ or $R = K_{R} - 0.8G$ by setting $\frac{\mathrm{d}R}{\mathrm{d}t} = 0$.

The G-nullclines are found to be $G = 0$ or $G = K_{G} - 0.09R$ by setting $\frac{\mathrm{d}G}{\mathrm{d}t} = 0$.

The equilibrium points are found from the intersections of the R and G nullclines, and are given to be: $(0,0)$, $(K_{R},0)$, $(K_{G},0)$ and $(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928})$.

The Jacobian matrix of the system is given by:

$$ J(R,G) = \begin{bmatrix} 0.61 - \frac{1.22}{K_{R}}R - \frac{0.488}{K_{R}}G & -\frac{0.488}{K_{R}} \\
-\frac{0.0738}{K_{G}}G & 0.82 - \frac{0.18}{K_{G}}G - \frac{0.0738}{K_{G}}R \end{bmatrix} $$

Now looking at the Jacobian at the equilibria points:

At the equilibrium point $(0,0)$, we have the Jacobian matrix

$$ J(0,0) = \begin{bmatrix} 0.61 & 0 \\
0 & 0.82 \end{bmatrix} $$

Since $J(0,0)$ is a diagonal matrix, it has two eigenvalues which are $\lambda_{1} = 0.61$ and $\lambda_{2} = 0.82$. Since these eigenvalues are non-positive real numbers, the corresponding fixed point is an unstable source;

At the equilibrium point $(K_{R},0)$, we have the Jacobian matrix

$$ J(K_{R},0) = \begin{bmatrix} -0.61 & -0.488 \\
0 & 0.82 - \frac{0.073K_{R}}{K_{G}}\end{bmatrix} $$

Since $J(K_{R},0)$ is an upper triangular matrix, it has two eigenvalues which are $\lambda_{1} = -0.61$ and $\lambda_{2} = 0.82 - \frac{0.073K_{R}}{K_{G}}$. Since these eignevalues are real numbers with $\lambda_{1} < 0 < \lambda_{2}$, the corresponding fixed point is an unstable saddlepoint.

At the equilibrium point $(0,K_{G})$, we have the Jacobian matrix

$$ J(0,K_{G}) = \begin{bmatrix} 0.61 - \frac{0.488K_{G}}{K_{R}} & 0 \\
-0.0738 & -0.82 \end{bmatrix} $$

Since $J(0,K_{G})$ is a lower triangular matrix, it has two eigenvalues which are $\lambda_{1} = 0.61 - \frac{0.488K_{G}}{K_{R}}$ and $\lambda_{2} = -0.82$. Since these eignevalues are real numbers with $\lambda_{2} < 0 < \lambda_{1}$, the corresponding fixed point is an unstable saddlepoint.

Finally at the equilibrium point $(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928})$, we have the Jacobian matrix

$$ J(\frac{K_{R}-0.8K_{G}}{0.928},\frac{K_{G}-0.09K_{R}}{0.928}) = \begin{bmatrix} \frac{0.488K_{G} - 0.61K_{R}}{0.928K_{R}} & \frac{0.3904K_{G} - 0.488K_{R}}{0.928K_{R}} \\
\frac{0.006642K_{R} - 0.0738K_{G}}{0.928K_{G}} & \frac{0.738K_{R} - 0.82K_{G}}{0.928K_{G}} \end{bmatrix} $$

$$
\lambda_{1,2} = \frac{1}{2} \left[
\frac{0.488 K_{G}^2 + 0.738 K_{R}^2 - 1.43 K_{R} K_{G}}{0.928 K_{R} K_{G}}
\pm
\sqrt{
\left(
\frac{0.488 K_{G}^2 + 0.738 K_{R}^2 - 1.43 K_{R} K_{G}}{0.928 K_{R} K_{G}}
\right)^2
- 4
\frac{
(0.488 K_{G} - 0.61 K_{R})(0.738 K_{R} - 0.82 K_{G})
- (0.3904 K_{G} - 0.488 K_R)(0.006642 K_{R} - 0.0738 K_{G})}
{(0.928)^2 K_{R} K_{G}}
}
\right]
$$

Since these eignevalues are real numbers with $\lambda_{2} < 0 < \lambda_{1}$, the corresponding fixed point is an unstable saddlepoint.


# Appendix

*Add any ideas/manuscripts, links and references below, treat this as a draft for now*

Jacobian formula

$$J(R,G) =\begin{bmatrix}
\frac{\partial{F}}{\partial{R}} & \frac{\partial{F}}{\partial{G}} \\
\frac{\partial{E}}{\partial{R}} & \frac{\partial{E}}{\partial{G}}
\end{bmatrix}$$

Where $F(R,G) = \frac{\mathrm{d}R}{\mathrm{d}t}$ and $E(R,G) = \frac{\mathrm{d}G}{\mathrm{d}t}$
