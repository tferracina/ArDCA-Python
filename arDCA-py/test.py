import numpy as np
import pytest
from numpy.testing import assert_allclose

from utils import sample_vectorized, entropy
from ar_types import ArNet

def generate_random_params(q, N):
    """Function to generate random parameters"""
    np.random.seed(42)

    # Generate p0
    p0 = np.random.rand(q)
    p0 /= p0.sum()

    # Generate H
    H = np.random.uniform(-1.0, 1.0, (N, q))

    # Generate J
    J = [np.random.uniform(-0.5,0.5, (q, q, i+1)) for i in range(N)]

    # Generate idx
    idxperm = np.random.permutation(N + 1)

    return H, J, p0, idxperm


def test_sample_function():
    np.random.seed(42)

    # define basic params
    q = 20    # number of states
    N = 5    # number of sites
    msamples = 10  # numer of sites to generate
    
    # Create an ArNet instance.
    H, J, p0, idxperm = generate_random_params(q, N)
    arnet = ArNet(H, J, p0, idxperm)

    # Run both sampling functions
    res_sample = sample_vectorized(arnet, msamples)

    # Print results.
    print("Output of sample (naïve loop-based):")
    print(res_sample)
    print("Shape:", res_sample.shape)

    # --- Basic tests ---
    expected_shape = (N + 1, msamples)
    assert res_sample.shape == expected_shape, "sample: Unexpected shape"
    # Check that all values are valid state indices (between 0 and q-1).
    assert np.all((res_sample >= 0) & (res_sample < q)), "sample: Values out of range"
    
    print("\nAll tests passed!")


def test_entropy():
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Define test parameters
    q = 3    # number of states
    N = 2    # number of sites
    msamples = 25  # number of samples
    

    Z = np.random.randint(0, q, size=(N, msamples))  
    W = np.ones(msamples) / msamples  # Equal weights
    
    # Test basic properties
    S = entropy(Z, W)
    assert len(S) == N
    assert np.all(S >= 0)
    assert np.all(S <= np.log(q))
    
    # Test with known distribution - adjusted to 0-based indexing
    Z_known = np.array([[0, 0, 1, 1, 2],
                        [1, 1, 0, 0, 2]])
    
    W_known = np.full(5, 0.2)  # Equal weights for 5 samples
    S_known = entropy(Z_known, W_known)
    
    # For uniform distribution over 2 states (p=0.4 each) and 1 state (p=0.2),
    # entropy should be: -0.4*log(0.4) - 0.4*log(0.4) - 0.2*log(0.2)
    expected_entropy = -0.4 * np.log(0.4) - 0.4 * np.log(0.4) - 0.2 * np.log(0.2)
    assert_allclose(S_known[0], expected_entropy, atol=1e-10)
    assert_allclose(S_known[1], expected_entropy, atol=1e-10)
    
    # Test with zero weights
    W_zero = np.zeros(msamples)
    W_zero[0] = 1.0
    S_zero = entropy(Z, W_zero)
    assert np.all(np.isclose(S_zero, 0.0))
    
    # Test with single state
    Z_single = np.zeros((N, msamples), dtype=int)  # Changed to use 0 instead of 1
    S_single = entropy(Z_single, W)
    assert np.all(np.isclose(S_single, 0.0))
    
    print("Entropy test results:")
    print("Random data entropy:", S)
    print("Known data entropy:", S_known)
    print("Single state entropy:", S_single)
    print("\nAll entropy tests passed!")

if __name__ == "__main__":
    #test_entropy()
    test_sample_function()
