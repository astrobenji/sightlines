import sightlines as los
import numpy as np

def test_halfway():
    z_mid, dist = los.find_midpoint(z_1=0, r_1=0, z_2=1, r_2=0)
    assert(z_mid == 0.5)
    assert(dist == 0.5)
    
def test_pythag():
    z_mid, dist = los.find_midpoint(z_1=0, r_1=4, z_2=6, r_2=4)
    assert(z_mid == 3)
    assert(dist == 5)
    
def test_no_soln():
    z_mid, dist = los.find_midpoint(z_1=2, r_1=4, z_2=2, r_2=40)
    assert(z_mid == None)
    assert(dist == None)
    
def test_triangles():
    # Exercise for the reader: draw this setup
    z_mid, dist = los.find_midpoint(z_1=0, r_1=1, z_2=3, r_2=2)
    assert(z_mid == 2.0)
    assert(dist  == np.sqrt(5))
    
def test_midpoint_too_short():
    z_mid, dist = los.find_midpoint(z_1=0, r_1=5, z_2=1, r_2=1)
    assert(z_mid == -11.5)
    