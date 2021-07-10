import sightlines as los
import numpy as np

BOX_SIZE = 100.0

def test_exclude_negative_z():
    origin = np.array([0,0,0])
    direction = np.array([1,0,0])
    gas_coords = np.array([[-1,4,0],[1,2,0]])
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(len(z_r_list) == 1)
    assert(z_r_list[0] == (1,2,1))
    
def test_keep_zero_z():
    origin = np.array([0,0,0])
    direction = np.array([1,0,0])
    gas_coords = np.array([[0,3,0],[2,3,4]])
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(len(z_r_list) == 2)
    assert(z_r_list[0] == (0,3,0))
    
def test_pythagoras():
    origin = np.array([0,0,0])
    direction = np.array([1,0,0])
    gas_coords = np.array([[0,3,0],[2,3,4]])
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(z_r_list[1] == (2,5,1))

def test_empty_list():
    origin = np.array([0,0,0])
    direction = np.array([1,0,0])
    gas_coords = []
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(z_r_list == [])
    
def test_origin_shift():
    origin = np.array([-2,-2,-2])
    direction = np.array([1,0,0])
    gas_coords = np.array([[0,3,-2],[2,1,2]])
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(len(z_r_list) == 2)
    assert(z_r_list[0] == (2,5,0))
    assert(z_r_list[1] == (4,5,1))

def test_direction_vector_too_big():
    origin = np.array([0,0,0])
    direction = np.array([4,0,0])
    gas_coords = np.array([[0,3,0],[2,3,4]])
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(len(z_r_list) == 2)
    assert(z_r_list[0] == (0,3,0))
    assert(z_r_list[1] == (2,5,1))
    
def test_direction_vector_too_small():
    origin = np.array([0,0,0])
    direction = np.array([0.5,0,0])
    gas_coords = np.array([[0,3,0],[2,3,4]])
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(len(z_r_list) == 2)
    assert(z_r_list[0] == (0,3,0))
    assert(z_r_list[1] == (2,5,1))

def test_periodic_boundary_corrections():
    origin = np.array([0,0,0])
    direction = np.array([1,0,0])
    gas_coords = np.array([[51,4,0],[49,2,0],[1,103,0],[-51,3,4]])
    z_r_list = los.get_z_r_list(origin, direction, gas_coords, BOX_SIZE)
    assert(len(z_r_list) == 3)
    assert(z_r_list[0] == (49,2,1))
    assert(z_r_list[1] == (1,3,2))
    assert(z_r_list[2] == (49,5,3))