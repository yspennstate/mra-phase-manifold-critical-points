# The moment-constrained alignment loss has critical points beyond the sign flips

Yitzchak Shmalo

In moment-constrained alignment for multireference alignment, introduced by Shahverdi, Ström and Andén
(arXiv:2409.04868), the loss is studied on the signals with the power spectrum of the true signal. They show that the
sign flips of the true signal are critical points and conjecture (Conjecture 3.8) that for signals of odd length these
are the only critical points, with the negative of the true signal as the global
maximizer. This paper gives a signal of length seven for which, at every noise level up to 1/100, another
signal on the same phase manifold has larger loss than every sign flip. The loss then has at least ten critical points up to
cyclic shifts, where the conjecture predicts eight, and by rescaling and a perturbation argument (Corollary 5) the conjecture fails at every positive
noise level on an open set of signals of length seven.

## Files

- `main.tex`, `main.pdf`: the paper.
- `check_counterexample.py`, `check_counterexample.out`: checks in exact rational arithmetic the intervals in the
  proof of Theorem 2 and the integer example in Remark 6. Standard library only.
- `monte_carlo.py`, `monte_carlo.out`: the Monte Carlo averages of Remark 7
  (NumPy, seed 1, 200000 noise samples).

## Running the checks

    python check_counterexample.py      # ends with ALL CHECKS PASS
    python monte_carlo.py               # a few seconds

## License

CC BY 4.0.
