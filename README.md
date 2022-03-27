# WishartMoments Package
## _A package to compute all invariant Wishart moments_


WishartMoments is a python Package using the mathematical software SAGE for computing all the invariant moments of a Wishart distributed variable.

## Features

- Compute the symbolic expression for the expectation of order k of of the variable W in terms of its parameters n and \Sigma.
- Evaluate the moments for concrete values of k, n, and \Sigma.
- Compute (symbolically) the moments for the inverse W^{-1} when possible, in terms of n, \Sigma and the size of the matrix.
- When usign the package inside a Jupyter Notebook, you can print the results as MathJax formulas to get a nice visualisation.
- Obtain the LaTex code that represents any of these expressions.
- You can compute the moments using the [Wishart Moment Calculator][WM] web interface.
 
## Installation
You should have already installed SAGE. If not, you can get it from [here][sage].

Open a SAGE shell and install WishartMoments using the SAGE pip utility by typing

```sh
sage -pip install WishartMoments
```
Once the installation is complete, to use the package in the command line you have to type `sage` in the SAGE shell to initiate SAGE, and once you see the prompt `sage:`, import the package (as any other Python package):
```sh
sage: import WishartMoments as wm
```

## How to use it

We will now show  how to compute the following Wishart moments of order $ 3 $:

\begin{equation*}
\begin{array}{lcl}
     \mathbb{E}(W^{3}) &=&  \left({n}^{3} + 3 \, {n}^{2} + 4 \, {n}\right) \Sigma^{3} + \left(2 \, {\left({n}^{2} + {n}\right)} {(\mathrm{tr} \, {\Sigma})}\right) \Sigma^{2} + \left({n} {(\mathrm{tr} \, {\Sigma})}^{2} + {\left({n}^{2} + {n}\right)} {(\mathrm{tr} \, {\Sigma}^{2})}\right) \Sigma \\
\mathbb{E}(W {(\mathrm{tr} \, W)}^{2}) &=& 8 \, {n} \Sigma^{3} + 4 \, {n}^{2} {(\mathrm{tr} \, {\Sigma})} \Sigma^{2} + \left({n}^{3} {(\mathrm{tr} \, {\Sigma})}^{2} + 2 \, {n}^{2} {(\mathrm{tr} \, {\Sigma}^{2})}\right) \Sigma
\end{array}
\end{equation*}
and

\begin{align*}
\mathbb{E}({{W^{-3}}})  &= \frac{{\left({n} - r - 1\right)} {{\Sigma}^{-3}}}{{\left({n} - r + 1\right)} {\left({n} - r\right)} {\left({n} - r - 3\right)} {\left({n} - r - 5\right)}} \frac{2 \, {{\Sigma}^{-2}} {(\mathrm{tr} \, {\Sigma}^{-1})}}{{\left({n} - r + 1\right)} {\left({n} - r\right)} {\left({n} - r - 3\right)} {\left({n} - r - 5\right)}} \\
&\phantom{=}+ \frac{{\left(2 \, {(\mathrm{tr} \, {\Sigma}^{-1})}^{2} + {n} {(\mathrm{tr} \, {\Sigma}^{-2})} - r {(\mathrm{tr} \, {\Sigma}^{-2})} - {(\mathrm{tr} \, {\Sigma}^{-2})}\right)} {{\Sigma}^{-1}}}{{\left({n} - r + 1\right)} {\left({n} - r\right)} {\left({n} - r - 1\right)} {\left({n} - r - 3\right)} {\left({n} - r - 5\right)}}
\end{align*}

After importing the Wishart moments package, we have to create an instance of the class `Expectations`.
```sh
sage: k=3
sage: expec = wm.Expectations(k)
```

We need to know how to reference the expressions of which we want to compute their expectations. We can get a list the expressions of order $ k $ by using the method `expressions` of `Expectations`, which returns a list of 2-element lists with the index of the portrait and the the expression for expectation of the moment corresponding to it.
```sh
sage: expec.expressions()
```
```
[0, W*tr(W, 1)^2]
[1, 2/3*W^2*tr(W, 1) + 1/3*W*tr(W, 2)]
[2, W^3]
```

Here `tr(A,j)` represents $ \mathrm{tr}(A^j) $. Therefore, to get $ W {(\mathrm{tr}\,{W})}^{2} $ we call the method `moment` with the index `0`

```
```

```sh
sage: expec.moment(0)
```

```
{
	'var': W*tr(W, 1)^2 ,
	'moment': 8*n*S^3 + 4*n^2*tr(S, 1)*S^2 + (n^3*tr(S, 1)^2 
	+ 2*n^2*tr(S, 2))*S
}
```

Similarly we use the index `2` to get $ W^3 $.
```sh
sage: expec.moment(2)
```
```
{
	'var': W^3 ,
	'moment': (n^3 + 3*n^2 + 4*n)*S^3 + (2*(n^2 + n)*tr(S, 1))*S^2 
	+ (n*tr(S, 1)^2 + (n^2 + n)*tr(S, 2))*S
}
```

As for the moment of the inverse, we can call `moment` with the argument `inverse` set to `True`.
```sh
sage: expec.expressions(inverse = True)
```
```
[0, inv(W, 1)*tr(W, -1)^2]
[1, 2/3*inv(W, 2)*tr(W, -1) + 1/3*inv(W, 1)*tr(W, -2)]
[2, inv(W, 3)]
```

Here `inv(A,j)` represents $ A^{-j} $. We use the index `2` to get $ W^{-3} $.
```sh
sage: expec.moment(2, inverse = True)
```
```
{
	'var': inv(W, 3) ,
	'moment': (n - r - 1)*inv(S, 3)/((n - r + 1)*(n - r)
	  *(n - r - 3)*(n - r - 5)) 
	+ 2*inv(S, 2)*tr(S, -1)/((n - r + 1)*(n - r)*(n - r - 3)*(n - r - 5)) 
	+ (2*tr(S, -1)^2 + n*tr(S, -2) - r*tr(S, -2) 
	- tr(S, -2))*inv(S, 1)/((n - r + 1)*(n - r)*(n - r - 1)
	  *(n - r - 3)*(n - r - 5))
}
```
### Latex code of the results
We can obtain the string with the $\LaTeX$ code representing these expressions by using the built-in function `latex`. For instance, if we want to get the code for the variable $W {(\mathrm{tr} \, W)}^{2}$, we should use the following commands
```sh
sage: latex(expec.moment(0)['var'])
```
```
W {(\mathrm{tr} \, W)}^{2}
```
and for its expectation, ${\Sigma} {n}^{3} {(\mathrm{tr} \, {\Sigma})}^{2} + 8 \, {\Sigma}^{3} {n}+ 2 \, {\left(2 \, {\Sigma}^{2} {(\mathrm{tr} \, {\Sigma})}+ {\Sigma} {(\mathrm{tr} \, {\Sigma}^{2})}\right)} {n}^{2}$,

```sh
sage: latex(expec.moment(0)['moment'])
```
```
8 \, {n} \Sigma^{3} + 4 \, {n}^{2} {(\mathrm{tr} \, {\Sigma})} \Sigma^{2}
+ \left({n}^{3} {(\mathrm{tr} \, {\Sigma})}^{2} 
+ 2 \, {n}^{2} {(\mathrm{tr} \, {\Sigma}^{2})}\right) \Sigma
```

Notice that an instance of the form `wm.Expectations(3)` only permits to compute moments of order $3$. To compute a moment of a different order, say $k=4$, the user has to instantiate a new object of the class `wm.Expectations(4)`. We will continue with the examples using the same parameter as we were doing so far, that is `k=3`, so that can keep using the same object `expec`.

### Evaluating the moments for specific values of the parameters

Now we show how to compute the numerical value of the moment $E(W \hbox{tr} (W^2)) $  for a Wishart distributions with parameters $n=10$ and 

\begin{equation*}
\Sigma = \begin{bmatrix}
4 & 1\\
1 & 3
\end{bmatrix}
\end{equation*}

>We want to remark that neither the package nor the website will check if the matrix $\Sigma$ is positive definite.

We first set the matrix $\Sigma$:

```sh
sage: Sigma = np.array([[4,1],[1,3]]);
```

To evaluate the moment we use the `evaluate_moment` where the parameters are `t` (the index of the expression in the list  `expec.expressions()` for which we require its expectation),  `n_param` (the numerical value for the parameter $n$), `Sigma` (the numerical value for the matrix $\Sigma$) and the boolean parameter `inverse` (`False` will compute the moment of $W$ and `True` the moment for $W^{-1}$). Here we need to compute $E(W \mathrm{tr} (W^2))$ and therefore
the parameters are passed as follows.

```sh
sage: ev = expec.evaluate_moment(t=0, n_param=10, Sigma=Sigma, inverse=False);
```

As `ev` is a dictonary, we can retrieve the variable by using the key `'var'`
```sh
sage: ev['var']
```
```
W*tr(W, 1)^2
```

To get the moment we use the key `'moment'`
 
```sh
sage: ev['moment']
```
```
array([[813600., 231120.],
       [231120., 582480.]], dtype=object)
```

## License

GNU Public License v3

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
   [sage]: <https://www.sagemath.org/>
   [WM]: <https://antunescarles.github.io/wishart-moments-calculator/>
