# sampler

[An ad hoc rejection sampler](https://github.com/YihaoGu/sampler/blob/main/rejection_sampler/rejection_sampler.R) is implemented for the target distribution
$$\pi(\eta) \propto \frac{\eta}{1+\eta^2} \exp \left(-a^2 \left(\eta-c\right)^2 \right),$$
where $a > 0, c \in \mathbb{R}$.

The modified ziggurat algorithm can be found [here](https://github.com/YihaoGu/sampler/blob/main/rejection_transformed/rejection_sampler_transformed.R). This implementation leverages the transformed target distribution
$$\pi(x) \propto \exp \left(-a^2 \left( \sqrt{e^{2x}-1} - c\right)^2\right),$$
where $x = \frac{1}{2}\log \left(1+\eta^2 \right)$.