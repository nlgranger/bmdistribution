# bmdistribution

An implementation of the multivariate Binary (Bernoulli) mixture model for 
Matlab. Can be used to model distribution of pixels in distributions of 
binary images for example.

## Usage

This library follows the matlab distribution class as closely as possible, and
more precisely the [Gaussian mixture
model](mathworks.com/help/stats/gmdistribution-class.html) one.

## Example

```matlab
Xapp = [bsxfun(@lt, rand(40000, 9), [.1 .2 .3 .4 .5 .6 .7 .8 .9]);
bsxfun(@lt, rand(10000, 9), [.9 .8 .7 .6 .5 .4 .3 .2 .1])];
bmm = fitbmdist(Xapp, 2, 'Options', struct('MaxIter', 1000, 'TolFun', 1e-4))
bmm.ComponentProportion
bmm.Means

Xval = [bsxfun(@lt, rand(1000, 9), [.1 .2 .3 .4 .5 .6 .7 .8 .9]);
        bsxfun(@lt, rand(1000, 9), [.9 .8 .7 .6 .5 .4 .3 .2 .1])];
clus = bmm.cluster(Xval);
# label numbers might be inverted, check numbers and verify confusion matrix.
```

## Limitations

Some methods are not implemented because they seemed less important.
Also the code has not been extensively tested yet
Please, contact me if you feel something is missing, or the code is not 
working as it should.

I can also merge pull requests if you wish to submit a fix or a new 
functionality.

## References

[JUAN, Alfons et VIDAL, Enrique. Bernoulli Mixture Models for Binary Images.
In : ICPR (3). 2004. p.
367-370.](http://users.dsic.upv.es/~ajuan/research/2004/Juan04_08b.pdf)
