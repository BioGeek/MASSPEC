"""The modern reconstruction gives the same answer in every process.

Python randomises string hashing per process, so code that iterates a set of
sequences can rank ties differently from one run to the next. The answer must
not depend on it.
"""
import os
import subprocess
import sys

PROGRAM = """
from masspec.evaluate import assemble
from masspec.masses import MASSES
from masspec.modern import FastObserved, run_extend
from masspec.simulate import simulate

UBIQUITIN = ("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDY"
             "NIQKESTLHLVLRLRGG")
simulated = simulate(UBIQUITIN, alkylated=False, seed=0)
result = run_extend(FastObserved(simulated.masses, 9.0), MASSES)
print(assemble(result.contigs))
"""


def _run_with_hash_seed(seed):
    environment = dict(os.environ, PYTHONHASHSEED=str(seed))
    completed = subprocess.run([sys.executable, "-c", PROGRAM],
                               capture_output=True, text=True, check=True,
                               env=environment)
    return completed.stdout


def test_same_answer_under_different_hash_seeds():
    assert _run_with_hash_seed(0) == _run_with_hash_seed(1)
