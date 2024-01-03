import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt

# My programme starts off by asking users for information concerning
# the number of cluster they want to define the scale then the
# programme try to create clusters.
# We put the point and abundance of them in a dictionary(dict_of_point),
# and produce a dictionary in which abundance of elements
# normalized for every points in the list. I also make a COPY of it.
# Then, we capture points from EVERY REMAINING POINTS in the COPY within
# the given "scale" defined and put them in a list assigned to that point,
# and combine both and put them in a dictionary meanwhile DELETE those
# captured points from the copy. Afterwards, for each points in the list,
# we try to append the adjacent list to them into the list which those points
# are in, by this way I can form lists of aggregated adjacent points and form
# clusters or outliers (where in their list only their own numbers like 1: [1])
# to form list of points in clusters. Finally, we split outliers from clusters
# by identifying the length of list, and find the centres of clusters by
# calculating means of coordinates, and by finding the length of cluster lists
# we can know how many points inside each clusters.
# after completing the main algorithm we enter a commanding centre where we can
# produce LIST OF CLUSTER GENERATED, DISTRIBUTION OF DISTANCE OF POINTS IN THE
# CLUSTER TO CENTRE, AVERAGE DISTANCE OF TWO CLUSTERS and AN 7-DIMENSIONS
# GRAPH OF OVERALL ELEMENT DISTIRBUTION


# this is a headline generater. The idea comes from our class
# will be used numerous time in the programme
def header():
    header = ""
    for i in range(60):
        header += "-"
    return header


# this is the start of my programme. It is a framework of my programme.
# first, we define the scale by either naming a number of cluster wished to get
# or defining a scale manually.
def scale_decider():
    while True:
        choice = input("press 'c' if you want to define a number of clusters\n"
                       "press 's' if you want to define scale yourself\n"
                       "input 'stop' if you want to stop\n")
        print(header())
        # if you want a number of cluster...
        if choice == "c":
            try:
                # you have to choose between 1 cluster to 16 clusters
                # the range of cluster is limited due to it works based on
                # some default scales I worked before. If the number of cluster
                # goes beyond 16, the clusters formed are very likely to be bad
                # e.g lots of clusters with less than 10 points inside
                scale_choice = abs(int(input("from 1 to 16, "
                                             "number of cluster you want?\n")))
                # if the input is not a integer
            except ValueError:
                print("check if your input is an integer")
                print(header())
            else:
                try:
                    scale = scale_default[scale_choice - 1]
                # if the number nomimated is out of range of 1 to 16
                except IndexError:
                    print("check your number of scale, index out of range")
                    print(header())
                else:
                    # otherwise, the input is legal and we start the programme
                    print(header())
                    # this is the main algorithm of data processing
                    data_generation(scale)
                    # this is the commanding centre of other features
                    function_choices()
                    # when the programme ends, it breaks the loop
                    # of this def function, then everything ends
                    break
        # if you want to set a scale manually...
        elif choice == "s":
            try:
                scale = float(input("please enter your scale\n"
                                    "a recommend scale is higher than 10,"
                                    " where there are little outliers\n"))
            # if your input is not a number
            except ValueError:
                print("check if your input is a number")
                print(header())
            else:
                # if your input is 0 or a negative number
                if scale <= 0:
                    print("this should be a positive number")
                else:
                    # start main programme, this is the same as above
                    print(header())
                    data_generation(scale)
                    function_choices()
                    break
        # if you want to stop beform starting programme
        elif choice == "stop":
            break
        # if the input to choice is not any of "c", "s" or "stop"
        else:
            print("input is the wrong")
            print(header())
    # at the end of programme we always draw a headline
    return print(header())


def data_generation(scale):
    # read file data
    original_file = open("cluster_data.csv", "r")
    data = original_file.readlines()
    # get rid of the headline
    data.pop(0)
    for keys in data:
        # we create a dictionary of point number and their info
        # using function called point_and_array()
        point_and_array(keys)
    for point in dict_of_points:
        # Normalize all the points' coordinates using function normalize()
        normalize(point)
    for point in dict_of_normalized_points:
        # now is the most important algorithm
        # we trap the adjacent points, details in below
        expansion(point, scale)
        # if the pool of points is empty, we end the process
        if len(dict_of_normalized_points_copy) != 0:
            continue
        else:
            break
    for start_point in dict_of_family:
        # another important process where we group the clusters
        cluster_collector(start_point)
    for cluster in cluster_dict:
        # we identfy the clusters from outliers
        # and calculate the centre of clusters
        cluster_mean(cluster)
    for outliers in list_of_outliers:
        # we pop the outliers from the cluster pool
        cluster_dict.pop(outliers)
    for centres_normal in list_of_centre_normal:
        # then we denormalize the centres' coordinates
        denormalize(centres_normal)
    for a in range(10):
        print("\n")
    # as we return, we end the main algorithm and get
    # into the commanding centre
    return


def point_and_array(keys):
    element_abundance_of_points = []
    keys = keys.strip()
    element_abundance_of_points_str = keys.split(",")
    # this is the name to the point, we define name as an int
    number = int(element_abundance_of_points_str.pop(0))
    # then we put elements abundance of the point inside a list
    for i in element_abundance_of_points_str:
        element_abundance_of_points.append(float(i))
    # and assign the list to the point name, and update it in
    # the dictionary of point
    dict_of_points.update({number: element_abundance_of_points})
    print("processing point and array for %s" % (number))
    # then, for the point we update the abundance of each elements
    # in the element_dict use this function
    element_dict_creation(element_abundance_of_points)
    return


def element_dict_creation(element_abundance_of_points):
    # the abundance of elements are added inside the lists of their
    # own element. the loop of adding based on each points so later when we
    # plot graph using matplotlib, they tells the info of scatters exactly as
    # the info of the same points
    for i in element_abundance_of_points:
        element = element_list[element_abundance_of_points.index(i)]
        element_dict[element].append(float(i))
    return


def normalize(point):
    print("processing normalize for %s" % (point))
    element_compile_in_point = dict_of_points[point]
    noramlized_abundance_list = []
    x = 0
    # normalization is done by min-max feature scaling
    for x in range(7):
        which_element = element_list[x]
        abundance_total = element_dict[which_element]
        dominator = float(element_compile_in_point[x]) - min(abundance_total)
        dinominator = max(abundance_total) - min(abundance_total)
        noramlized_result = (dominator / dinominator) * 100
        noramlized_abundance_list.append(noramlized_result)
        element_dict_normalized[element_list[x]].append(noramlized_result)
    # create a new array of abundances and assign them to the point name
    # the same as in point_and_array
    noramlized_abundance_array = np.array(noramlized_abundance_list)
    # and update the point and their arrays into dict_of_normalized_point
    # similar to dic_of_point
    dict_of_normalized_points.update({point: noramlized_abundance_array})
    # we also make a copy of it. we use it to delete adjacent points later.
    dict_of_normalized_points_copy.update({point: noramlized_abundance_array})
    return


def expansion(point, scale):
    # we try to use a 7-dimensional hypersphere to trap the adjacent points
    # we identify the "radius" to the hypersphere
    initial_scale = np.zeros(7)
    upper_scale = initial_scale + scale
    # we trap the points using point_spotter function
    # and finally put the original point with its adjacent point list
    # inside dic_of_family
    dict_of_family.update({point:
                           point_spotter(point, upper_scale)})
    return


def point_spotter(point, upper_scale):
    point_array = np.array(dict_of_normalized_points[point])
    dist_scale = distance.euclidean(upper_scale, np.zeros(7))
    new_adjacent_list = []
    removal_list = []
    for other_point in dict_of_normalized_points_copy:
        # identify adjacent points from every points left in the point pool
        coordinate = np.array(dict_of_normalized_points_copy[other_point])
        dist = abs(distance.euclidean(coordinate, point_array))
        # if we identify one
        if dist <= dist_scale:
            print("adjacent spotted around %s" % (point))
            new_adjacent_list.append(other_point)
            # we remove the spotted adjacent point from the pool
            # I do this because I want to keep every points appear
            # inside adjacent list of any points only once.
            # it is good for grouping in the families later
            # only outliers cannot be captured by any points so they could
            # be left in the pool and get captured only by themselves
            # so we won't make any misjudgements about adjacence
            removal_list.append(other_point)
    # we do the removing
    for i in removal_list:
        dict_of_normalized_points_copy.pop(i)
    # return the point list
    return new_adjacent_list


def cluster_collector(start_point):
    # the first cluster put in as a cluster
    if len(cluster_dict) == 0:
        cluster_dict.update({start_point: dict_of_family[start_point]})
    else:
        judge_list = []
        for point_exist in cluster_dict:
            for other_adjacent in cluster_dict[point_exist]:
                judge_list.append(bool(other_adjacent == start_point))
                # if the point is in any existing clusters
                if bool(other_adjacent == start_point) is True:
                    # we append the adjacent list to that point to the cluster
                    new_adjacent = dict_of_family[start_point]
                    cluster_dict[point_exist].extend(new_adjacent)
                    print("point %s added in a cluster" % (new_adjacent))
                else:
                    continue
        # if the point is not in any of the existeing clusters
        # it is either a new cluster or an outlier
        # but right now they are the similar thing to us
        # so it is added to the cluster dictionary
        # to clearify, outliers are like 1:[1], while points that look like
        # 1:[2] is actually an adjacent to some other clusters since
        # their own points have been captured, so they themselves should also
        # be captured above
        if judge_list.count(True) == 0:
            print("create new centre from %s" % (start_point))
            cluster_dict.update({start_point: dict_of_family[start_point]})
    return


def cluster_mean(cluster):
    sum_of_points = np.zeros(7)
    length = len(cluster_dict[cluster])
    print(length)
    if length == 1:
        # in this case, the points have length 1 can only be outliers
        list_of_outliers.append(cluster)
    else:
        # we calculate the centre of cluster using mean of
        # coordinates of all the points inside the cluster
        for point in cluster_dict[cluster]:
            point_value = dict_of_normalized_points[point]
            sum_of_points += point_value
        centre_new = sum_of_points / length
        # the new centre is transformed into arrays and append to a list
        list_of_centre_normal.append(np.array(centre_new))
    return


def denormalize(centres_normal):
    # this is a reverse process of normalization in above
    print("denormalizing centre")
    element_abundance_centre = []
    centres_normal = list(centres_normal)
    for i in centres_normal:
        index = centres_normal.index(i)
        element_overall_2 = element_dict[element_list[index]]
        dinominator_2 = max(element_overall_2) - min(element_overall_2)
        minimum_abundance = min(element_overall_2)
        denormaled_abuncance = (i / 100) * dinominator_2 + minimum_abundance
        element_abundance_centre.append(denormaled_abuncance)
    list_of_centre.append(np.array(element_abundance_centre))
    return


# this is the commanding centre
def function_choices():
    # we print the cluster info automatically at the start
    # at the same time it renames the clusters
    cluster_info()
    while True:
        print("please select the function you want to use")
        function = input("press 'a' for list of clusters generated\n"
                         "press 'b' for average distance between two cluster\n"
                         "press 'c' for distribution of point distance to"
                         " centre in cluster\n"
                         "press 'd' for distribution of 7 elements together\n"
                         "press 'stop' to stop the programme\n")
        if function == "a":
            cluster_info()
        elif function == "b":
            average_distance_finder()
        elif function == "c":
            distribution_inside_cluster()
        elif function == "d":
            plot_graph()
        elif function == "stop":
            break
        else:
            print("input error, try again")
            print(header())
    return


def cluster_info():
    # print out the abundace of elements in centres
    cluster_key = list(cluster_dict.keys())
    total_points = 0
    for blank in range(3):
        print("\n")
    print(header())
    print("******cluster centres and numbers in the cluster******")
    for key_of_cluster in cluster_key:
        # this is a calcualtion of total points collected inside all clusters
        # if my programme works correctlt the total number, added with
        # number of outliers, should be 1027, same as number of points
        # in csv file.
        # print_out is another function returns number of points, and itself
        # print the abundance of each elements of centre
        total_points += print_out(key_of_cluster, cluster_key)
    print(header())
    print("total number of cluster spotted is %g" % (len(list_of_centre)))
    print(header())
    print("total points collected in clusters is %g" % (total_points))
    print(header())
    # print the number of outliers and their point numbers
    if len(list_of_outliers) == 0:
        print("there is no outliers")
    else:
        print("there are %g outliers" % (len(list_of_outliers)))
        print(list_of_outliers)
    print(header())
    return print(header)


def print_out(key_of_cluster, cluster_key):
    index_in_cluster = cluster_key.index(key_of_cluster) + 1
    print(header())
    # we reassign the clusters with numerical order of name
    print("Centre %s" % (index_in_cluster))
    point_compile = cluster_dict[key_of_cluster]
    cluster_dict.pop(key_of_cluster)
    # now reassign name by modifying the cluster dict, put the new name to it
    cluster_dict.update({index_in_cluster: point_compile})
    point_around_cluster = len(point_compile)
    # print out abundance
    for x in range(7):
        print("abundance of %s: %.3f" %
              (element_list[x], list_of_centre[index_in_cluster - 1][x]))
    print("number of points in this cluster is %g" % (point_around_cluster))
    if point_around_cluster <= 5:
        print("Not a good cluster, try to use a"
              "bigger scale or a smaller cluster number")
    # return the number of points inside the cluster
    return point_around_cluster


# this is the function finding the distribution inside a given cluster
def distribution_inside_cluster():
    while True:
        print(header())
        print("distribution of distance inside cluster")
        # remind you the number of cluster captured so you won't mistype
        print("total number of clusters is %g" % (len(list_of_centre)))
        print(header())
        cluster_selection = input("please enter centre number or 'stop'\n")
        # in case you want to stop
        if cluster_selection == "stop":
            print(header())
            break
        else:
            try:
                cluster_selection = int(cluster_selection)
            # if your input is not a number
            except ValueError:
                print("value error, try again")
            else:
                try:
                    # create histogram with function histogram_creation
                    histogram_creation(cluster_selection)
                    print("histogram of cluster centre %g is created" %
                          (cluster_selection))
                # if your input is not within the total number of clusters
                except IndexError:
                    print("check your input, your input should be less"
                          "than or equal to the total number of clusters")


def histogram_creation(cluster_selection):
    # create histogram with raw data, no further data processes
    # if the clusters have distributions similar to
    # a normal distribution or a bimordal pattern, it is a good cluster.
    # in fact, if the selected cluster number is small or the selected
    # scale is relatively large, this pattern should be more obvious
    distance_list = []
    selected_centre = list_of_centre[cluster_selection - 1]
    selected_cluster_list = cluster_dict[cluster_selection]
    for point in selected_cluster_list:
        point_array = np.array(dict_of_points[point])
        distance_within = distance.euclidean(selected_centre, point_array)
        distance_list.append(distance_within)
    distance_array = np.array(distance_list)
    plt.hist(distance_array)
    plt.title("distribution of distance to center in cluster %g"
              % (cluster_selection))
    plt.xlabel("distance to centre")
    plt.ylabel("number of points in ranges")
    plt.grid()
    return plt.show()


def average_distance_finder():
    print(header())
    print("average distance finder")
    print("total number of clusters is %g" % (len(list_of_centre)))
    print(header())
    while True:
        # when naming both centres, you can input stop and it will ends
        centre_1 = input("the number of first centre or 'stop'\n")
        print(header())
        centre_2 = input("the number of second centre or 'stop'\n"
                         "if you have already selected stop above,"
                         " press enter\n")
        if centre_1 == "stop" or centre_2 == "stop":
            print(header())
            break
        else:
            try:
                centre_1 = int(centre_1)
                centre_2 = int(centre_2)
            # in case you put a wrong data type
            except ValueError:
                print("value error try again")
                print(header())
            else:
                # if your centre number is out of range
                if centre_1 > len(list_of_centre)\
                   or centre_2 > len(list_of_centre):
                    print("number centre out of range, try again")
                    print(header())
                else:
                    average = average_distance(centre_1, centre_2)
                    print(header())
                    print("the average distance between centre %g"
                          " and centre %g is %.4f" %
                          (centre_1, centre_2, average[0]))
                    print(header())
                    print("discrepancy in abundance")
                    for abundance in average[1]:
                        new_index = list(average[1]).index(abundance)
                        print("abundance in %s: %.4f"
                              % (element_list[new_index], abundance))
                    print(header())
    return


def average_distance(centre_1, centre_2):
    # for distance I put the euclidean distance
    # for every points in one cluster, I calculate the distance of
    # it to every point in the other cluster, sum up the distances and
    # calculate the average
    # it does the same for abundance of individual elements
    list_for_distance = []
    list_for_chemistry = []
    for points in cluster_dict[centre_1]:
        coordinate = dict_of_points[points]
        for points_compare in cluster_dict[centre_2]:
            coordinate_comp = dict_of_points[points_compare]
            # calculate distance use euclidean distance
            their_distance = distance.euclidean(coordinate, coordinate_comp)
            list_for_distance.append(their_distance)
            # calculate difference in abundance
            abundance_difference = np.array(np.array(coordinate)
                                            - np.array(coordinate_comp))
            list_for_chemistry.append(abundance_difference)
    # for the distance
    average_distance_clusters = sum(list_for_distance) / len(list_for_distance)
    # for difference in abundance
    abundance_difference = np.array(sum(list_for_chemistry)
                                    / len(list_for_chemistry))
    return [average_distance_clusters, abundance_difference]


def plot_graph():
    # this is the scatter points that mixed up 7 dimension datas together
    # the interpretation is in below, it will show up if you choose to
    # plot the graph.
    print(header)
    print("loading 7 dimension graph, read the instruction below")
    print(header())
    print("here is a instruction:\n"
          "x, y and z axis represent abundance of Fe, Al and Mg\n"
          "size of point (small to large), represent abundance "
          "of Mn (small to large)\n"
          "colour of point (red to purple in regular spectrum order"
          "and finally black), represent abundance of Si (small to large)\n"
          "transparency (transparent to opague),"
          " represent abundance of Ni (small to large)\n"
          "shape (trangle_left, trangle_right,"
          " square, pentagon, hexagon1, diamond, x) in sequence,"
          " represent abundance of Co (small to large")
    print(header())
    colour_list = ["r", "y", "g", "c", "b", "m", "k"]
    shape_list = ["<", ">", "s", "p", "h", "D", "X"]
    ax = plt.axes(projection="3d")
    element = element_dict.keys()
    for i in element:
        element_dict[i] = np.array(element_dict[i])
    for x in range(1027):
        colour = element_dict_normalized["Si"][x] / 100
        colour_from_Si = colour_list[int((colour * 7) - 1)]
        alpha = element_dict_normalized["Ni"][x] / 100
        shape_sequence = element_dict_normalized["Co"][x] / 100
        shape_from_Co = shape_list[int((shape_sequence * 7) - 1)]
        ax.scatter3D(element_dict["Fe"][x], element_dict["Al"][x],
                     element_dict["Mg"][x], s=100*element_dict["Mn"][x],
                     c=colour_from_Si, alpha=alpha, marker=shape_from_Co)
    ax.set_xlabel("abundance of Fe")
    ax.set_ylabel("abundance of Al")
    ax.set_zlabel("abundance of Mg")
    return plt.show()


# default scale for different number of clusters, max is 16 clusters
# if cluster sets beyond 16, the clusters we get may not be representative
scale_default = [45, 41, 35, 33, 25, 22, 15, 13.5, 13,
                 12.3, 12.2, 12, 11.93, 11.9, 11.5, 11]
# dictionaries about points
dict_of_points = {}
dict_of_normalized_points = {}
# this dictionary is solely for checking if every
dict_of_normalized_points_copy = {}
# dictionaries of "family trees"
dict_of_family = {}
# dictionary of collecting abundance for plot
element_dict = {"Fe": [], "Al": [], "Mg": [],
                "Mn": [], "Si": [], "Ni": [], "Co": []}
element_dict_normalized = {"Fe": [], "Al": [], "Mg": [],
                           "Mn": [], "Si": [], "Ni": [], "Co": []}
# dictionary of collecting cluster
cluster_dict = {}
# list of compiled centre arrays, will be used in point and array
list_of_centre_normal = []
list_of_centre = []
# list of outlier points
list_of_outliers = []
# element list in sequence, will be used many times
element_list = list(element_dict.keys())
# colour and shape file, will be used in plot_graph
colour_list = ["r", "y", "g", "c", "b", "m", "k"]
shape_list = ["<", ">", "s", "p", "h", "D", "X"]
# now we start the main programme, see details in the front
scale_decider()
