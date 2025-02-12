import numpy as np
import pytest
from numpy.testing import assert_allclose

from utils import sample, entropy
from ar_types import ArNet

def test_sampling_functions():
    np.random.seed(42)

    # define basic params
    q = 3    # number of states
    N = 2    # number of sites
    msamples = 5  # numer of sites to generate

    p0 = np.array([0.2, 0.5, 0.3])   # initial prob dist over states

    # H list of len N where each entry is a vector of len q
    H = np.array([[ 0.1, 0.2, 0.3],
                  [-0.1, 0.0, 0.1]])

    # J: list of len N
    # For site 0, J[0] should have shape (q, q, 1)
    J0 = np.array([
        [[0.01], [0.02], [0.03]],
        [[0.04], [0.05], [0.06]],
        [[0.07], [0.08], [0.09]]
    ])
    # For site 1, J[1] should have shape (q, q, 2)
    J1 = np.array([
        [[0.10, 0.11], [0.12, 0.13], [0.14, 0.15]],
        [[0.16, 0.17], [0.18, 0.19], [0.20, 0.21]],
        [[0.22, 0.23], [0.24, 0.25], [0.26, 0.27]]
    ])
    J = [J0, J1]

    # idxperm: a permutation of length (N+1). This will be used to reorder the rows.
    idxperm = np.array([2, 0, 1])
    
    # Create an ArNet instance.
    arnet = ArNet(H, J, p0, idxperm)

    # Run both sampling functions
    res_sample = sample(arnet, msamples)
    #res_vectorized = sample_vectorized(arnet, msamples)

    # Print results.
    print("Output of sample (naïve loop-based):")
    print(res_sample)
    print("Shape:", res_sample.shape)
    
    #print("\nOutput of sample_vectorized:")
    #print(res_vectorized)
    #print("Shape:", res_vectorized.shape)

    # --- Basic tests ---
    expected_shape = (N + 1, msamples)
    assert res_sample.shape == expected_shape, "sample: Unexpected shape"
    #assert res_vectorized.shape == expected_shape, "sample_vectorized: Unexpected shape"
    
    # Check that all values are valid state indices (between 0 and q-1).
    assert np.all((res_sample >= 0) & (res_sample < q)), "sample: Values out of range"
    #assert np.all((res_vectorized >= 0) & (res_vectorized < q)), "sample_vectorized: Values out of range"
    
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
    test_entropy()
