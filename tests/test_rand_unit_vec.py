import sightlines as los
import numpy as np

def test_vec_length():
    vec = los.generate_rand_unit_vec()
    assert(len(vec) == 3)

def test_vec_norm():
    vec = los.generate_rand_unit_vec()
    np.testing.assert_almost_equal(np.linalg.norm(vec), 1.0)