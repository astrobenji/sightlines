'''
A package designed to work with Illustris/TNG, to find the absorption 
metallicity seen in a galaxy from a GRB/QSO along a line of sight.
'''

__author__ = ['Benjamin Metha', ]
__date__ = '2019-12-17'
__cite__ = 'https://github.com/astrobenji'

import numpy as np
import random

###########################################################
#                                                         #
#           The main two (3) functions that human         #
#                   beings should use                     #
#                                                         #
###########################################################
    
def get_los_metallicity_and_H1_fraction(origin, direction, gas_dict, ionising_distance=0.0, box_size = 75000.0):
    # Taking advantage of the cylindrical symmetry of this problem, reduce it 
    # to 2D.
    z_r_list = get_z_r_list(origin, direction, gas_dict['Coordinates'], box_size, ionising_distance)
    # Find out how long the line of sight spends in each of our gas cells.
    segments = compute_len_in_each_cell(z_r_list)
    
    metal_density = 0
    col_density   = 0
    H1_density    = 0
    for cell_id in segments:
        metal_density += segments[cell_id]*gas_dict['Density'][cell_id]*gas_dict['GFM_Metallicity'][cell_id]
        col_density += segments[cell_id]*gas_dict['Density'][cell_id]
        H1_density += segments[cell_id]*gas_dict['Density'][cell_id]*gas_dict['NeutralHydrogenAbundance'][cell_id]*(1-gas_dict['GFM_Metallicity'][cell_id])
    
    # Avoiding divide by zero error.
    if col_density == 0:
        return 0
    
    los_Z = metal_density/col_density
    return (los_Z, H1_fraction)

def get_los_metallicity(origin, direction, gas_dict, ionising_distance=0.0, box_size = 75000.0):
    '''
    The only real function anyone outside of this package needs to use.
    Gets the metallicity of a galaxy with gas cells given by gas_dict, as it
    would be seen by a GRB originating at origin and travelling in direction.
    
    Parameters
    ----------
    origin: (3,) numpy array
        The location of the cell from which the GRB is produced
        
    direction: (3,) numpy array
        The unit vector of the direction of the GRB
        
    gas_dict: dict-type.
        Contains all the details about our gas cells.
        MUST CONTAIN AT LEAST:
         * `Coordinates`, which gives the coordinates of the gas cells within 
           our periodic box, in units of ckpc/h;
         * `Density`, which gives the density of gas in this cell in units of
           (10^10M_sun/h)/(ckpc/h)^3;
         * `GFM_Metallicity` - The ratio M_Z/M_tot for this cell, where M_Z is
           the total mass all metal elements (above He)
        Names and units come directly from https://www.tng-project.org/data/docs/specifications/
        Nothing bad will happen if the dict has other fields - this program
        will not edit the gas_dict.
    
    ionising_distance: float
        The distance from the origin of the GRB to the first part of the LOS for which 
        cells are ionised. Should be set to 0.1 **physical** kpc (Vreeswijk+12). Kept at 
        0 so you don't get lazy and just assume it's at the right level already.
       
    box_size: float
        The size of our periodic box. Defaulted to 75000 - this is the value 
        in ckpc/h for TNG100.
                
    Returns
    -------
    los_Z: float
        The metallicity seen along this line of sight.
    
    '''
    # Taking advantage of the cylindrical symmetry of this problem, reduce it 
    # to 2D.
    z_r_list = get_z_r_list(origin, direction, gas_dict['Coordinates'], box_size, ionising_distance)
    # Find out how long the line of sight spends in each of our gas cells.
    segments = compute_len_in_each_cell(z_r_list)
    
    metal_density = 0
    col_density = 0
    for cell_id in segments:
        metal_density += segments[cell_id]*gas_dict['Density'][cell_id]*gas_dict['GFM_Metallicity'][cell_id]
        col_density += segments[cell_id]*gas_dict['Density'][cell_id]
    
    # Avoiding divide by zero error.
    if col_density == 0:
        return 0
    
    los_Z = metal_density/col_density
    return los_Z
    
def mc_get_los_metallicity(origin, gas_dict, n_trials=100, ionising_distance=0.0, box_size=75000.0):
    '''
    Same as get_los_metallicity, but it repeats it `n_trials` times with
    different directions. Ideal for seeing how angular variation can change
    the observed metallicity from the same GRB host star.
    
    Parameters
    ----------
    origin: (3,) numpy array
        The location of the cell from which the GRB is produced
        
    gas_dict: dict-type.
        Contains all the details about our gas cells.
        MUST CONTAIN AT LEAST:
         * `Coordinates`, which gives the coordinates of the gas cells within 
           our periodic box, in units of ckpc/h;
         * `Density`, which gives the density of gas in this cell in units of
           (10^10M_sun/h)/(ckpc/h)^3;
         * `GFM_Metallicity` - The ratio M_Z/M_tot for this cell, where M_Z is
           the total mass all metal elements (above He)
        Names and units come directly from https://www.tng-project.org/data/docs/specifications/
        Nothing bad will happen if the dict has other fields - this program
        will not edit the gas_dict.
        
    n_trials: int
        The number of times to repeatedly generate lines of sight from this
        origin point.
    
    ionising_distance: float
        The distance from the origin of the GRB to the first part of the LOS for which 
        cells are ionised. Should be set to 0.1 **physical** kpc (Vreeswijk+12). Kept at 
        0 so you don't get lazy and just assume it's at the right level already.
    
    box_size: float
        The size of our periodic box. Defaulted to 75000 - this is the value 
        in ckpc/h for TNG100.
        
    Returns
    -------
    los_Z: (n_trials,) ndarray
        The metallicity seen along `n_trials` lines of sight that all start
        from the origin cell.
    
    '''
    los_Z = np.zeros(n_trials)
    for ii in range(n_trials):
        direction = generate_rand_unit_vec()
        los_Z[ii] = get_los_metallicity(origin, direction, gas_dict, box_size, ionising_distance)
    return los_Z
    
###########################################################
#                                                         #
#                      Subprocesses                       #
#                  (in order of appearance)               #
#                                                         #
###########################################################

def get_z_r_list(origin, direction, gas_coords, box_size, ionising_distance=0.0):
    '''
    Create a list of gas cells with their distance from the LOS and their 
    location along the LOS.
    
    Parameters
    ----------
    origin: (3,) numpy array
        The location of the cell from which the GRB is produced
        
    direction: (3,) numpy array
        The unit vector of the direction of the GRB
    
    gas_coords: (3,N) numpy array
        Contains the location in 3-space of each gas cell in a periodic box.
    
    box_size: float
        The size of our periodic box.  
    
    ionising_distance: float
        The distance from the origin of the GRB to the first part of the LOS for which 
        cells are ionised. Should be set to 0.1 **physical** kpc (Vreeswijk+12). Kept at 
        0 so you don't get lazy and just assume it's at the right level already.  
    
    Returns
    -------
    z_r_list: list of 3-tuples.
        These tuples will contain:
        z - the position of each gas cell along the LOS;
        r - the distance of this gas cell from the LOS; and
        cell_id - the index of the gas cell in question in this sheet.
                  Allows the cell's other properties to be looked up later.
    '''
    # Check that direction is a unit vector -- if not, fix this.
    if np.linalg.norm(direction) != 1.0:
        direction = direction/np.linalg.norm(direction)
        
    if len(gas_coords) == 0:
        return []
    
    rel_coords = gas_coords - origin
    
    # ** Removing periodic box effects **
    # If a coord is over 0.5*box_size, shift everything back box_size
    # If a coord is less than -0.5*box_size, move everything over box_size
    # This will ensure all coords lie between -0.5*box_size and 0.5*box_size, 
    # and nothing straddles an edge.
    rel_coords = rel_coords + (rel_coords < -0.5*box_size)*box_size - (rel_coords > 0.5*box_size)*box_size
    
    z_r_list = []
    for cell_id, cell_posn in enumerate(rel_coords):
        # Find z-distance by taking the dot product (scalar projection) of 
        # the cell's position relative to the origin and the direction of the LOS.
        z = np.dot(direction, cell_posn)
        # If a cell originates from behind the burst/within the ionising 
        # distance, do not include it.
        if z<ionising_distance:
            continue
        r = np.linalg.norm(cell_posn - z*direction)
        z_r_list.append( (z,r, cell_id) )
    
    return z_r_list

def compute_len_in_each_cell(z_r_list):
    '''
    Computes the length of the sightline through each cell. Saves output
    as a dict.
    
    Parameters
    ----------
    z_r_list: list of 3-tuples
        These tuples will contain:
        z - the position of each gas cell along the LOS;
        r - the distance of this gas cell from the LOS; and
        cell_id - the index of the gas cell in question in this sheet.
                  Allows the cell's other properties to be looked up later.
    
    Returns
    -------
    segment_dict: dict
        Keys are the IDs of each cell through which this line of sight passes.
        Values are the length of the line segment through that cell.
    
    '''
    if z_r_list == []:
        return {}
    
    # Sort the z_r_list by z!
    z_r_list.sort()
    
    # Put first cell and last cell into a queue.
    n_cells = len(z_r_list)
    z_start = 0 # For QSOs this should be z_r_list[0][0]
    z_end   = z_r_list[n_cells-1][0]
    start_cell, _ = closest_cell_to(z_start, z_r_list)
    end_cell, _ = closest_cell_to(z_end, z_r_list)
    queue = [(z_start, z_end, start_cell, end_cell)]
    
    # Container for tuples consisting of (0) the ID of each cell and (1) the
    # length of the segment of the sightline that stays within this gas cell.
    segment_dict = {}
    
    while queue:
        z_start, z_end, start_cell, end_cell = queue.pop(0)
        
        z_start_cell, r_start_cell, start_cell_id = z_r_list[start_cell]
        z_end_cell, r_end_cell, end_cell_id       = z_r_list[end_cell]
        
        # If the closest point to the start of this segment==the closest point
        # to the end of the segment, then the entire segment will lie within
        # the same Voronoi cell. This comes from the fact that Voronoi cells
        # are convex.
        if start_cell == end_cell:
            # Add this segment to our dict
            if start_cell_id not in segment_dict.keys():
                segment_dict[start_cell_id] = 0
            segment_dict[start_cell_id] += z_end - z_start
            continue
        
        # Find the point on the sightline that is equally distant from these 
        # two cells. Save this, together with this distance
        midpoint, dist_to_midpoint = find_midpoint(z_1=z_start_cell, 
                                                   r_1=r_start_cell, 
                                                   z_2=z_end_cell, r_2=r_end_cell)
        
        if midpoint == None:
            print("Error: The midpoint does not exist.")
            exit(1)
        elif midpoint < z_start:
            print("Error: Midpoint too early.")
            exit(1)
        elif midpoint > z_end:
            print("Error: Midpoint too late.")
            exit(1)
        
        sublist = z_r_list[start_cell+1:end_cell]
        closest_cell_sublist_index, closest_cell_dist = closest_cell_to(midpoint, sublist)
        
        if closest_cell_sublist_index == None: 
            # Then there is no cells between these two.
            # Add these two segments to our segment dict...
            if start_cell_id not in segment_dict.keys():
                segment_dict[start_cell_id] = 0
            segment_dict[start_cell_id] += midpoint - z_start
            
            if end_cell_id not in segment_dict.keys():
                segment_dict[end_cell_id] = 0
            segment_dict[end_cell_id] += z_end - midpoint
            
            # ...and move on to the next item in the queue.
            continue
        else:
            # Convert this sublist index into a full list index.
            closest_cell = closest_cell_sublist_index + start_cell + 1

        if closest_cell_dist < dist_to_midpoint:
            # Then split the line into two smaller lines, and repeat this 
            # process for each.
            queue.append((z_start, midpoint, start_cell, closest_cell))
            queue.append((midpoint, z_end, closest_cell, end_cell))
        else:
            # Then there are no other cells closer to the midpoint.
            # Add these two segments to our dict.
            if start_cell_id not in segment_dict.keys():
                segment_dict[start_cell_id] = 0
            segment_dict[start_cell_id] += midpoint - z_start
            
            if end_cell_id not in segment_dict.keys():
                segment_dict[end_cell_id] = 0
            segment_dict[end_cell_id] += z_end - midpoint
            
    return segment_dict

def find_midpoint(z_1, r_1, z_2, r_2):
    '''
    Finds the z-value of the point on the line r=0 that is equally distant 
    between the cell at (z_1, r_1) and (z_2, r_2).
    
    Parameters
    ----------
    z_1: float
        z-coord of the first cell.
    r_1: float
        r-coord of the first cell.
    z_2: float
        z-coord of the second cell.
    r_2: float
        r-coord of the second cell.
        
    Returns
    -------
    z_mid: float
        z-coord of the midpoint between these cells on the line r=0.
        
    dist: float
        Distance of this point to either cell.
    '''
    if z_1 == z_2:
        # Then there is no midpoint between these cells.
        return None, None
    
    # I used the quadratic formula here so you don't have to :~)
    z_mid = 1.0*(r_2**2 + z_2**2 - r_1**2 - z_1**2)/(2*(z_2 - z_1))
    # Then a bit of Pythagoras...
    dist  = np.sqrt(r_1**2 + (z_1 - z_mid)**2)
    
    return z_mid, dist

def closest_cell_to(point, z_r_list):
    '''
    Finds the closest cell to any given point on the line, given our z_r_list.
    
    Parameters
    ----------
    point: float
        The place on the line we are trying to find the closest point to.
    
    z_r_list: list of 3-tuples.
        These tuples will contain:
        z - the position of each gas cell along the LOS;
        r - the distance of this gas cell from the LOS; and
        cell_id - the index of the gas cell in question in this sheet.
                  Allows the cell's other properties to be looked up later.
    
    Returns
    -------
    closest_cell_index: int
        The ID of the cell that is closest to our given point.
    
    closest_cell_dist: float
        The distance of this closest cell to the specified point.
    
    '''
    if z_r_list==[]:
        return None, None
    
    z_vals = np.array([tup[0] for tup in z_r_list])
    r_vals = np.array([tup[1] for tup in z_r_list])
    
    d2_vals = r_vals**2 + (z_vals - point)**2
    
    closest_cell_index = np.argmin(d2_vals)
    closest_cell_dist  = np.sqrt(d2_vals[closest_cell_index])
    
    return closest_cell_index, closest_cell_dist

def generate_rand_unit_vec():
    '''
    Generates random points on a 2-Sphere. Uses algorithm described in https:/
    /math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-
    vector-in-3d
    
    Returns
    -------
    vec: (3,) numpy array
    
    A randomly pointing unit vector (Cartesian Coordinates).
    '''
    z     = random.uniform(-1,1)
    theta = random.uniform(0, 2*np.pi)
    
    rho    = np.sqrt(1 - z**2)
    vec = np.array([rho*np.cos(theta), rho*np.sin(theta), z])
    
    return vec 
    
def rotate_around_axis(v,erot,angle):
    rotmeasure=np.linalg.norm(erot)
    erot=erot/rotmeasure;
    norme=np.dot(v,erot)
    vplane=v-norme*erot
    plnorm=np.linalg.norm(vplane)
    ep=vplane/plnorm
    eo=np.cross(erot,ep)
    vrot=(np.cos(angle)*ep+np.sin(angle)*eo)*plnorm+norme*erot
    return(vrot)
    