import sightlines as los
import numpy as np

def test_halfway():
    short_z_r_list = [(0,0,1), (1,0,10)]
    seg_dict = los.compute_len_in_each_cell(short_z_r_list)
    assert(len(seg_dict)==2)
    assert(seg_dict[1]==0.5)
    assert(seg_dict[10]==0.5)

def test_equal_10():
    nice_z_r_list = [(0,1,10), (1,1,11), (2,1,12), (3,1,13), (4,1,14), 
                     (5,1,15), (6,1,16), (7,1,17), (8,1,18), (9,1,19)]
    seg_dict = los.compute_len_in_each_cell(nice_z_r_list)
    assert(len(seg_dict)==10)
    assert(seg_dict[10] == 0.5)
    assert(seg_dict[11] == 1.0)
    assert(seg_dict[12] == 1.0)
    assert(seg_dict[13] == 1.0)
    assert(seg_dict[14] == 1.0)
    assert(seg_dict[15] == 1.0)
    assert(seg_dict[16] == 1.0)
    assert(seg_dict[17] == 1.0)
    assert(seg_dict[18] == 1.0)
    assert(seg_dict[19] == 0.5)

def test_middle_tie():
    pythagoras_hill = [(0,0,10), (9,3,21), (15,3,22), (24,0,88)]
    seg_dict = los.compute_len_in_each_cell(pythagoras_hill)
    assert(len(seg_dict)==4)
    assert(seg_dict[10] == 5)
    assert(seg_dict[21] == 7)
    assert(seg_dict[22] == 7)
    assert(seg_dict[88] == 5)
    
def test_4_way_tie():
    circular_hill = [(0,0,10), (2,1,21), (3,2,22), (5,2,24), (6,1,25), (8,0,76)]
    seg_dict = los.compute_len_in_each_cell(circular_hill)
    assert(len(seg_dict)==4)
    assert(seg_dict[10]==1.25)
    assert(seg_dict[21]==2.75)
    assert(seg_dict[25]==2.75)
    assert(seg_dict[76]==1.25)
    
def test_worst_case_for_efficiency():
    meano_zeno = [(1,0,1), (2,0,2), (4,0,4), (8,0,8), (16,0,16), (32,0,32), (64,0,64)]
    seg_dict = los.compute_len_in_each_cell(meano_zeno)
    assert(len(seg_dict)==7)
    assert(seg_dict[1]==1.5)
    assert(seg_dict[2]==1.5)
    assert(seg_dict[4]==3)
    assert(seg_dict[8]==6)
    assert(seg_dict[16]==12)
    assert(seg_dict[32]==24)
    assert(seg_dict[64]==16)
    
def test_no_galaxy():
    just_the_origin = [(0,0,0)]
    seg_dict = los.compute_len_in_each_cell(just_the_origin)
    assert(len(seg_dict)==1)
    assert(seg_dict[0]==0)
    
def test_empty():
    seg_dict = los.compute_len_in_each_cell([])
    assert(len(seg_dict)==0)

def test_the_thing_I_was_worried_about_before_but_will_probably_be_totally_fine_since_I_redesigned_my_code():
    calm_down = [(0,0,0), (9,5,4), (9,3,1)]
    seg_dict = los.compute_len_in_each_cell(calm_down)
    assert(len(seg_dict)==2)
    assert(seg_dict[0]==5)
    assert(seg_dict[1]==4)
    
def test_ends_are_not_closest():
    high_ends = [(0,100,100), (1,1,10), (9,1,20), (10,100, 200)]
    seg_dict = los.compute_len_in_each_cell(high_ends)
    assert(len(seg_dict)==2)
    assert(seg_dict[10]==5)
    assert(seg_dict[20]==5)
    