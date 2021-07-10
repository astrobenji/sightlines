import sightlines as los
import numpy as np

def test_closest_int():
    nice_z_r_list = [(0,0,10), (1,0,11), (2,0,12), (3,0,13), (4,0,14), 
                     (5,0,15), (6,0,16), (7,0,17), (8,0,18), (9,0,19)]
    cell_index, dist = los.closest_cell_to(3.1, nice_z_r_list)
    assert(cell_index == 3)
    np.testing.assert_almost_equal(dist, 0.1)
    
def test_tie_breaker():
    nice_z_r_list = [(0,0,10), (1,0,11), (2,0,12), (3,0,13), (4,0,14), 
                     (5,0,15), (6,0,16), (7,0,17), (8,0,18), (9,0,19)]
    cell_index, dist = los.closest_cell_to(3.5, nice_z_r_list)
    assert(cell_index == 3) # In the event of a tie, the first min is chosen.
    np.testing.assert_almost_equal(dist, 0.5)
    
def test_empty():
    cell_index, dist = los.closest_cell_to(3.5, [])
    assert(cell_index == None)
    assert(dist == None)