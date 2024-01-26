import numpy as np 
from matplotlib import pyplot as plt
from scipy.spatial.distance import euclidean

###### 
def Calculate_Distances(coordinate_tuple):
    N = len(coordinate_tuple)
    distances = []
    for i in range(N-1):
        x0, y0, z0 = coordinate_tuple[i][1], coordinate_tuple[i][2], coordinate_tuple[i][3]
        x1, y1, z1 = coordinate_tuple[i+1][1], coordinate_tuple[i+1][2], coordinate_tuple[i+1][3]
        d = np.sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)
        distances.append(d)
    return distances

def Calculate_Angular_Distances(coordinate_tuple):
    N = len(coordinate_tuple)
    angular_distances = []
    for i in range(N-2):
        l1, m1, n1 = coordinate_tuple[i][1]-coordinate_tuple[i+1][1], coordinate_tuple[i][2]-coordinate_tuple[i+1][2], coordinate_tuple[i][3]-coordinate_tuple[i+1][3]
        l2, m2, n2 = coordinate_tuple[i+1][1]-coordinate_tuple[i+2][1], coordinate_tuple[i+1][2]-coordinate_tuple[i+2][2], coordinate_tuple[i+1][3]-coordinate_tuple[i+2][3]
        cos_t = np.abs(l1*l2+m1*m2+n1*n2) / ((np.sqrt(l1*l1+m1*m1+n1*n1))*(np.sqrt(l2*l2+m2*m2+n2*n2)))
        angular_distances.append(cos_t)
    return angular_distances

def Calculate_Distances_All(coordinate_tuple):
    N = min([len(coordinate_tuple), 200]) #limit to first 100 grains
    distances = []
    for i in range(N):
        x0, y0, z0 = coordinate_tuple[i][1], coordinate_tuple[i][2], coordinate_tuple[i][3]
        for j in range(N):
            if (i!=j):
                x1, y1, z1 = coordinate_tuple[j][1], coordinate_tuple[j][2], coordinate_tuple[j][3]
                d = np.sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)
                distances.append(d)
    return distances

def Calculate_Angular_Distances_All(coordinate_tuple):
    N = min(len(coordinate_tuple), 25)
    angular_distances = []
    for i in range(N):
        for j in range(N):
            if (j!=i):
                l1, m1, n1 = coordinate_tuple[i][1]-coordinate_tuple[j][1], coordinate_tuple[i][2]-coordinate_tuple[j][2], coordinate_tuple[i][3]-coordinate_tuple[j][3]
                for k in range(N):
                    if (k!=j and k!=i):
                        l2, m2, n2 = coordinate_tuple[k][1]-coordinate_tuple[j][1], coordinate_tuple[k][2]-coordinate_tuple[j][2], coordinate_tuple[k][3]-coordinate_tuple[j][3]
                        cos_t = np.abs(l1*l2+m1*m2+n1*n2) / ((np.sqrt(l1*l1+m1*m1+n1*n1))*(np.sqrt(l2*l2+m2*m2+n2*n2)))
                        angular_distances.append(cos_t)
    return angular_distances


######

def Check_Impact_Parameter_Match(coordinate_tuple1, theta1, phi1, coordinate_tuple2):
    N1 = len(coordinate_tuple1)
    # 1 -> 2
    x0, y0, z0 = coordinate_tuple1[N1-1][1], coordinate_tuple1[N1-1][2], coordinate_tuple1[N1-1][3]
    x1, y1, z1 = coordinate_tuple2[0][1], coordinate_tuple2[0][2], coordinate_tuple2[0][3]
    theta0, phi0 = theta1, phi1
    dz = z0 - z1 
    x0 = x0 + dz*np.tan(theta0)*np.cos(phi0)
    y0 = y0 + dz*np.tan(theta0)*np.sin(phi0)
    dx, dy = x0-x1, y0-y1
    b_for = np.sqrt(dx*dx + dy*dy)
    return b_for


def Check_Distance_Match(coordinate_tuple1, coordinate_tuple2):
    N1 = len(coordinate_tuple1)
    N2 = len(coordinate_tuple2)
    # 1 -> 2
    x0, y0, z0 = coordinate_tuple1[N1-1][1], coordinate_tuple1[N1-1][2], coordinate_tuple1[N1-1][3]
    x1, y1, z1 = coordinate_tuple2[0][1], coordinate_tuple2[0][2], coordinate_tuple2[0][3]
    d = np.sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)

    l1, m1, n1 = coordinate_tuple1[N1-1][1]-coordinate_tuple1[0][1], coordinate_tuple1[N1-1][2]-coordinate_tuple1[0][2], coordinate_tuple1[N1-1][3]-coordinate_tuple1[0][3]
    l2, m2, n2 = coordinate_tuple2[N2-1][1]-coordinate_tuple2[0][1], coordinate_tuple2[N2-1][2]-coordinate_tuple2[0][2], coordinate_tuple2[N2-1][3]-coordinate_tuple2[0][3]
    cos_t = np.abs(l1*l2+m1*m2+n1*n2) / ((np.sqrt(l1*l1+m1*m1+n1*n1))*(np.sqrt(l2*l2+m2*m2+n2*n2)))

    return d, cos_t


def Convert_to_TX_TY(theta, phi):
    tan_theta, tan_phi = np.tan(theta), np.tan(phi)
    tx = np.sqrt((tan_theta*tan_theta)/(1+tan_phi*tan_phi))
    if (phi>=np.pi/2 and phi<=1.5*np.pi):
        tx = -tx
    ty = np.sqrt(tan_theta*tan_theta - tx*tx)
    if (phi>=np.pi):
        ty=-ty
    return tx, ty

def Convert_to_Theta_Phi(thetaX, thetaY):
    theta = np.arctan(thetaX**2 + thetaY**2)
    phi = np.arctan(thetaY/thetaX)
    return theta, 2*np.pi-phi


def merge_and_remove_duplicates(tracks1, tracks2):
    merged_tracks = tracks1 + tracks2

    # Create a set to keep track of seen tracks and a result list for tracks3
    seen_tracks = set()
    tracks3 = []

    for track in merged_tracks:
        if track not in seen_tracks:
            tracks3.append(track)
            seen_tracks.add(track)

    return tracks3


def extract_tracks_from_couples(couple):

    N = len(couple)
    tuple_size = 2 
    # Split the original tuple into N tuples
    tracks_list = [couple[i:i + tuple_size] for i in range(0, N, tuple_size)]
    return tracks_list


def check_for_duplicates(list1, list2):
    set1 = set(list1)
    set2 = set(list2)

    # Use set intersection to find common elements
    common_elements = set1.intersection(set2)

    if common_elements:
        return True  # There are duplicates
    else:
        return False  # No duplicates found
    

def find_duplicates(list1, list2):
# Create sets from the tuples in both lists
    set1 = set(list1)
    set2 = set(list2)

    # Find the common elements (duplicates) between the two sets
    duplicates = set1.intersection(set2)

    return list(duplicates)


def eliminate_equivalent_lists(lists):
    # Create a dictionary to store the longest list for each unique tuple
    max_lists = {}

    for lst in lists:
        # Extract the unique tuple from the list
        unique_tuple = set([item[0] for item in lst])
        
        # Convert the tuple to a frozenset to make it hashable
        unique_tuple = frozenset(unique_tuple)
        
        # Check if this tuple is already in the dictionary
        if unique_tuple in max_lists:
            # If the current list is longer than the one in the dictionary, replace it
            if len(lst) > len(max_lists[unique_tuple]):
                max_lists[unique_tuple] = lst
        else:
            # If the tuple is not in the dictionary, add it with the current list
            max_lists[unique_tuple] = lst

    # Convert the values of the dictionary back to a list
    result = list(max_lists.values())

    return result


def sort_tuples_by_z(tuples_list):
    sorted_tuples = sorted(tuples_list, key=lambda x: x[2], reverse=True)
    return sorted_tuples


def line_fit_and_properties(points, PLOTS=0):
    points_array = np.array(points)

    # Perform a linear fit for the YZ plane
    coefficients_yz, _, _, _, _ = np.polyfit(points_array[:, 1], points_array[:, 2], 1, full=True)

    # Extract slope for the YZ plane
    slope_yz = coefficients_yz[0]
    intercept_yz = coefficients_yz[1]

    # Calculate the angle in the YZ plane (arctan(slope_yz))
    angle_yz = np.arctan(slope_yz)

    # Perform a linear fit for the XZ plane
    coefficients_xz, _, _, _, _ = np.polyfit(points_array[:, 0], points_array[:, 2], 1, full=True)

    # Extract slope for the XZ plane
    slope_xz = coefficients_xz[0]
    intercept_xz = coefficients_xz[1]

    # Calculate the angle in the XZ plane (arctan(slope_xz))
    angle_xz = np.arctan(slope_xz)

    # Plotting
    if (PLOTS):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Plotting the points
        ax1.scatter(points_array[:, 1], points_array[:, 2], label='Points')
        ax1.set_xlabel('Y')
        ax1.set_ylabel('Z')
        ax1.set_title('Points in YZ plane')

        # Plotting the YZ plane fit line
        fit_line_yz = np.polyval([slope_yz, intercept_yz], points_array[:, 1])
        ax1.plot(points_array[:, 1], fit_line_yz, color='red', label='Fit Line')
        ax1.legend()

        # Plotting the points
        ax2.scatter(points_array[:, 0], points_array[:, 2], label='Points')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Z')
        ax2.set_title('Points in XZ plane')

        # Plotting the XZ plane fit line
        fit_line_xz = np.polyval([slope_xz, intercept_xz], points_array[:, 0])
        ax2.plot(points_array[:, 0], fit_line_xz, color='red', label='Fit Line')
        ax2.legend()

        plt.show()

    # Estimate the length of the fit line between the first and last points in 3D space
    length_3d = euclidean(points_array[0, :], points_array[-1, :])

    return angle_yz, angle_xz, length_3d


def FindGrains(total_grains_list, view, mtid, close_view_list):
    for i in range(len(total_grains_list[close_view_list.index(view)])):
        if (total_grains_list[close_view_list.index(view)][i][0] == mtid):
            return total_grains_list[close_view_list.index(view)][i][1]
        
def GetGrainsIntervals(grains_list):
    xlist, ylist, zlist = [], [], []
    for entry in grains_list:
        xlist.append(entry[0])
        ylist.append(entry[1])
        zlist.append(entry[2])
    return xlist, ylist, zlist


def TestIntervals(int1, int2):
    min1, min2 = min(int1), min(int2)
    max1, max2 = max(int1), max(int2)

    if ((min1 < min2 and max1 > max2) or (min2 < min1 and max2 > max1)):
        return True 
    else:
        return False
    


def CheckImpactN(s1, s2, dzMax=6000.):
    imp = 1e+10
    p1 = np.array([s1.X() - dzMax * s1.TX(), s1.Y() - dzMax * s1.TY(), s1.Z() - dzMax])
    p2 = np.array([s1.X() + dzMax * s1.TX(), s1.Y() + dzMax * s1.TY(), s1.Z() + dzMax])
    p3 = np.array([s2.X() - dzMax * s2.TX(), s2.Y() - dzMax * s2.TY(), s2.Z() - dzMax])
    p4 = np.array([s2.X() + dzMax * s2.TX(), s2.Y() + dzMax * s2.TY(), s2.Z() + dzMax])

    pa = np.zeros(3)
    pb = np.zeros(3)
    mua, mub = 0.0, 0.0

    if LineLineIntersect(p1, p2, p3, p4, pa, pb, mua, mub):
        parallel = False
        pv = 0.5 * (pa + pb)
        imp = np.linalg.norm(pa - pb)
    else:
        parallel = True
        pv = np.array([0.5 * (s1.X() + s2.X()), 0.5 * (s1.Y() + s2.Y()), 0.5 * (s1.Z() + s2.Z())])
        inside = 0
        imp = 2. * DistancePointLine3(pv, p1, p2, inside)

    return imp, pv, parallel


def LineLineIntersect(p1, p2, p3, p4, pa, pb, mua, mub):
    p13 = np.subtract(p1, p3)
    p43 = np.subtract(p4, p3)

    EPS = 1.E-6

    if np.all(np.abs(p43) < EPS):
        return False

    p21 = np.subtract(p2, p1)

    if np.all(np.abs(p21) < EPS):
        return False

    d1343 = np.dot(p13, p43)
    d4321 = np.dot(p43, p21)
    d1321 = np.dot(p13, p21)
    d4343 = np.dot(p43, p43)
    d2121 = np.dot(p21, p21)

    denom = d2121 * d4343 - d4321 * d4321

    if np.abs(denom) < EPS:
        return False

    numer = d1343 * d4321 - d1321 * d4343

    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    pa[0] = p1[0] + mua * p21[0]
    pa[1] = p1[1] + mua * p21[1]
    pa[2] = p1[2] + mua * p21[2]
    pb[0] = p3[0] + mub * p43[0]
    pb[1] = p3[1] + mub * p43[1]
    pb[2] = p3[2] + mub * p43[2]

    return True

def DistancePointLine3(p, l1, l2, inside):
    l1l2 = l2 - l1
    l1p = p - l1

    t = np.dot(l1p, l1l2) / np.dot(l1l2, l1l2)

    if t < 0.0:
        t = 0.0
        inside = False
    elif t > 1.0:
        t = 1.0
        inside = False
    else:
        inside = True

    projection = l1 + t * l1l2
    distance = np.linalg.norm(p - projection)

    return distance


def CheckProb(s1, s2, eDZmax, eZbin, eAbin, zpos1, zpos2):

    c = 1
    dz = abs(s1.Z() - s2.Z())
    if dz > eDZmax:
        c = 0

    isign = 0
    dtx, dty, deltaZ = 0, 0, 0

    if zpos1 != zpos2:
        deltaZ = dz + eZbin
        isign = -1
    else:
        deltaZ = eDZmax - dz / 2.0
        isign = 1

    dtx = abs(s1.TX() - isign * s2.TX()) + eAbin
    if abs(s1.X() - s2.X()) > dtx * deltaZ:
        c = 0

    dty = abs(s1.TY() - isign * s2.TY()) + eAbin
    if abs(s1.Y() - s2.Y()) > dty * deltaZ:
        c = 0
    return c
