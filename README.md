# sampler

[rejection_sampler](./rejection_sampler) implements an ad hoc sampler for the target distribution
$$\pi(\eta) \propto \frac{\eta}{1+\eta^2} \exp \left(-a^2 \left(\eta-c\right)^2 \right),$$
where $a > 0, c \in \mathbb{R}$.

[rejection_transformed](./rejection_transformed/) implements the modified ziggurat algorithm based on the transformed target distribution
$$\pi(x) \propto \exp \left(-a^2 \left( \sqrt{e^{2x}-1} - c\right)^2\right),$$
where $x = \frac{1}{2}\log \left(1+\eta \right)^2$.