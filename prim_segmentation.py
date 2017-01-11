from PIL import Image
import numpy as np
from scipy.ndimage import filters
import random
import time

# Converting a colour image into gray scale image and pre-processing step for 
# initializing dictionaries for Union find function.
def segment_prim(image, k_val, sigma, min_size):
    start_time = time.time()
    im = np.array(Image.open(image).convert('L'))
    row, col = im.shape
    im = filters.gaussian_filter(im, sigma)
    threshold = {}
    vertices = {}
    another_v = {}
    num = 0
    for r_1 in range(row):
        for c_1 in range(col):
            vertices[(r_1, c_1)] = num
            another_v[num] = (r_1, c_1)
            threshold[num] = threshold_calc(1, k_val)
            num += 1
        threshold = {}
    for i in range(row * col):
        threshold[i] = threshold_calc(1, k_val)

    collect_vertices = list(range(0, row * col))
    u = Unionfind(collect_vertices)

    collect_edges = {}
    var = 0
    for c in range(1, col - 1):
        for r in range(1, row - 1):
            edge = neighbor_edge(r, c)
            for e in edge:
                d = diff(im, e[0], e[1], r, c)
                collect_edges[var] = [vertices[r, c], vertices[e[0], e[1]], d]
                a = u.find(vertices[(r, c)])
                b = u.find(vertices[(e[0], e[1])])
                if a != b:
                    if d <= threshold[a] and d <= threshold[b]:
                        u.union(a, b)
                        threshold[a] = d + threshold_calc(len(u.set_items[u.set_id[a]]), k_val)
    for i in collect_edges:
        a = u.find(collect_edges[i][0])
        b = u.find(collect_edges[i][1])
        if (a != b) and ((len(u.set_items[a]) < min_size) or (len(u.set_items[b]) < min_size)):
            u.union(a, b)

    colors = {}
    for i in u.set_items:
        colors[i] = (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))

    im = np.array(Image.open(image))
    for c in u.set_items:
        for v in u.set_items[c]:
            r_c, v_c = another_v[v]
            im[r_c, v_c] = colors[c]
    im = Image.fromarray(im)
    # saving the output image
    im.save("tomato-prim-segment.png")
    time_program = time.time() - start_time
    print("computing the segments", time_program)
    return time_program


# computing edges along the way
def neighbor_edge(row_val, col_val):
    down = (row_val-1, col_val)
    up = (row_val+1, col_val)
    left = (row_val, col_val-1)
    right = (row_val, col_val+1)
    return [up, down, left, right]

# region growing to group image pixels whose intensities are lesser than the threshold calculated using threshold function.
class Unionfind:
    def __init__(self, i_s):
        self.set_id = dict((s, s) for s in i_s)
        self.set_items = dict((s, [s]) for s in i_s)

    def find(self, element):
        return self.set_id[element]

    def union(self, s1, s2):
        set1 = self.set_items[s1]
        set2 = self.set_items[s2]
        if len(set1) >= len(set2):
            short_set = s2
            long_set = s1
        else:
            short_set = s1
            long_set = s2
        for item in self.set_items[short_set]:
            self.set_id[item] = long_set
            self.set_items[long_set].append(item)
        del self.set_items[short_set]

# calculating the edge weights while converting image into a graph.
def diff(image_val, w, h, new_w, new_h):
    im_weight = np.abs(int(image_val[w, h]) - int(image_val[new_w, new_h]))
    return im_weight

# threshold function pivotal in calculating segments 
def threshold_calc(size, k):
    return k/size


def main():
    sigma = 0.8
    k_val = 300
    image = 'tomatoes-grape.jpg'
    min_size = 500
    segment_prim(image, k_val, sigma, min_size)

if __name__ == "__main__":
    main()

