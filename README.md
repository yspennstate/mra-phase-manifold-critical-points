# The moment-constrained alignment loss has critical points beyond the sign flips

Yitzchak Shmalo

Shahverdi, Ström and Andén (arXiv:2409.04868) restrict the expected template-alignment loss of multireference
alignment to the signals with the power spectrum of the true signal, and conjecture (Conjecture 3.8) that for signals
of odd length the only critical points are the sign points, with the negative of the true signal as the global
maximum. This note gives a signal of length seven for which, at every noise standard deviation up to 1/100, another
signal on the same torus has larger loss than every sign point. The loss then has at least ten critical orbits under
cyclic shifts, where the conjecture predicts eight, and by rescaling the conjecture fails at every positive noise
level.

## Files

- `main.tex`, `main.pdf`: the paper.
- `check_counterexample.py`, `check_counterexample.out`: checks in exact rational arithmetic the intervals in the
  proof of Theorem 2 and the integer example in Remark 5. Standard library only.
- `monte_carlo.py`, `monte_carlo.out`: the Monte Carlo averages of Remark 6
  (NumPy, seed 1, 200000 noise samples).

## Running the checks

    python check_counterexample.py      # ends with ALL CHECKS PASS
    python monte_carlo.py               # a few seconds

## License

CC BY 4.0.
