import numpy as np

from utils import sample
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

test_sampling_functions()