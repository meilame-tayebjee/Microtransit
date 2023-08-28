import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import random
import sys

RADIUS = int(sys.argv[1])
SAMPLE_SIZE = int(sys.argv[2])
DELTA = int(sys.argv[3])
DELTA2 = int(sys.argv[4])
iter = int(sys.argv[5])


def simulatePoints(size, center, radius):
    """size : INT, center = 2-dim array (coordinates), radius = float"""
    rhos = radius *np.sqrt(np.random.random(size = size))
    thetas = np.random.uniform(0, 2*np.pi, size = size)
    xs = center[0] + rhos * np.cos(thetas)
    ys = center[1] + rhos * np.sin(thetas)

    return xs, ys


def list_pairs(n):
    pairs = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            if j != i:
                pairs.append((i, j))
    return pairs


def exhaustive_arcfile(locations):
    size = len(locations)
    pairs = list_pairs(size)
    begin = [i[0] for i in pairs]
    end = [i[1] for i in pairs]
    distances = []
    for pair in pairs:
            i, j = pair
            distances.append(np.sqrt(np.sum((locations.loc[i] - locations.loc[j])**2)))

    return begin, end, distances

def reshuffle_array(arr):
    res = np.copy(arr)
    for i in range(len(arr)):
        ind = int(np.random.choice([j for j in range(len(arr)) if j != i]))
        res[i] = arr[ind]
    return res


def simulate_setting(RADIUS, SAMPLE_SIZE, DELTA1, DELTA2):
    """
    This function simulates a setting with 4 stations and a central station
    Around each station, simulation of SAMPLE_SIZE origins and SAMPLE_SIZE destinations, uniformly distributed in a circle of radius RADIUS
    DELTA is the distance between the two stations on the horizontal axis, DELTA2 is the distance between the two stations on the vertical axis
    Each origin is assigned randomly to a destination
    """

    R = RADIUS
    SIZE = SAMPLE_SIZE
    # Simulate clients location uniformly around a center ("Station")
    STATION = [0,0]
    xs_init, ys_init = simulatePoints(2*SIZE // 3, STATION, R)
    # Simulate clients destination uniformly around another center ("CENTRAL")
    CENTRAL = [DELTA1, 0]
    xs_central, ys_central = simulatePoints(SIZE // 3, CENTRAL, R)
    xs_central_new, ys_central_new = simulatePoints(SIZE // 3, CENTRAL, R)

    TERMINUS = [DELTA1 + DELTA2, 0]
    xs_final, ys_final = simulatePoints(2*SIZE // 3, TERMINUS, R)

    xs_dest = np.concatenate((xs_central, xs_final))
    ys_dest = np.concatenate((ys_central, ys_final))
    xs_origin = np.concatenate((xs_init, xs_central_new))
    ys_origin = np.concatenate((ys_init, ys_central_new))

    filename = "Linear_Delta_" + str(DELTA) + "_Delta2_" + str(DELTA2) + "_Radius_" + str(RADIUS) + "_SampleSize_" + str(SAMPLE_SIZE) + "_iter_" + str(iter) +".png"

    #Plotting and saving figure
    plt.figure(figsize = (10,5))
    plt.scatter(xs_origin, ys_origin, color = 'b')
    plt.scatter(xs_central_new, ys_central_new, color = 'b')
    plt.scatter(xs_central, ys_central,  color = 'g')
    plt.scatter(xs_final, ys_final,  color = 'g')
    plt.scatter(STATION[0], STATION[1], color = 'r', label = 'Station')
    plt.plot([STATION[0], CENTRAL[0]], [STATION[1], CENTRAL[1]], color ='r')
    plt.plot([STATION[0], TERMINUS[0]], [STATION[1], TERMINUS[1]], color ='r',  label = 'Transit line')
    plt.scatter(CENTRAL[0], CENTRAL[1], color = 'r', label = 'Central Station')
    plt.scatter(TERMINUS[0], TERMINUS[1], color = 'r', label = 'Terminus')
    plt.plot([xs_origin, xs_dest], [ys_origin, ys_dest], color = 'r', linestyle = ":", alpha = 0.5)
    plt.legend()
    plt.savefig('visualizations/setting_visualizations/' + filename)

    centers = [STATION, CENTRAL, TERMINUS]
    return xs_origin, ys_origin, xs_dest, ys_dest, centers

def generate_OnD_csv(xs_origin, ys_origin, xs_dest, ys_dest, centers, SAMPLE_SIZE):

    xs = np.concatenate((xs_origin, xs_dest))
    ys = np.concatenate((ys_origin, ys_dest))
    location_df = pd.DataFrame(data = {"location" : range(1, 2*SAMPLE_SIZE+1), "lat_coord" : xs, "long_coord" : ys})
    location_df = location_df.set_index("location")

    for i, center in enumerate(centers):
        location_df.loc[len(location_df.index) + 1] = center
    
    begin, end, distances = exhaustive_arcfile(location_df)
    arcfile = pd.DataFrame({"startloc" : begin, "endloc":end, "traveltime":[d for d in distances], "distance":distances})
    arcfile = arcfile.set_index("startloc")

    location_df.to_csv("/Users/Meilame/Desktop/MIT/data/locations_" + str(SAMPLE_SIZE)+ ".csv")
    arcfile.to_csv("/Users/Meilame/Desktop/MIT/data/arcs_" + str(SAMPLE_SIZE)+ ".csv")

    return 

def generate_Transit_csv(centers):
    location_df = pd.DataFrame(data = {"location" : range(1, len(centers) + 1), "lat_coord" : [center[0] for center in centers], "long_coord" : [center[1] for center in centers]})
    location_df = location_df.set_index("location")
    location_df.to_csv("/Users/Meilame/Desktop/MIT/data/transit_locations.csv")
    begin, end, distances = exhaustive_arcfile(location_df)
    arcfile = pd.DataFrame({"startloc" : begin, "endloc":end, "traveltime":distances, "distance":distances})
    arcfile = arcfile.set_index("startloc").iloc[[0, 3]]
    arcfile.to_csv("/Users/Meilame/Desktop/MIT/data/transit_arcs.csv")

    return

def main():
    xs_origin, ys_origin, xs_dest, ys_dest, centers = simulate_setting(RADIUS, SAMPLE_SIZE, DELTA, DELTA2)

    generate_OnD_csv(xs_origin, ys_origin, xs_dest, ys_dest, centers, SAMPLE_SIZE)
    generate_Transit_csv(centers)
    print("Succesfully generated setting")
    return

if __name__ == "__main__":
    main()



