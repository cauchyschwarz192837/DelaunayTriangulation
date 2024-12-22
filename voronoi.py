from primitives import *
from delaunay import *
import random

def get_voronoi(T, margin=1.2):
    '''given a Triangulation T, get the set of vertices and sets of 
    bounded and semi-infinite edges of the dual of T, i.e., the Voronoi diagram
    of T's points, T.pts.'''
    
    tris = T.get_triangles()
    
    verts = set()
    bounded_edges = set()
    semi_infinite_edges = set()
    
    # for each pair of triangles that share an edge,
    #   draw their dual vertices and connect them with dual edge
    for ta, adj_tris in tris.items():
        ca = ta.circum().get_center()
        verts.add(ca)

        for tb in adj_tris:
        
            cb = tb.circum().get_center()    
            verts.add(cb)

            bounded_edges.add(Segment(ca,cb))

    # rest of this is a hack to draw semi-infinite edges
    #   as finite-length since matplotlib does not support them natively.
    
    # goal is to "clip" the semi-infinite edges within a bounding box whose
    #   size is determined by the given "margin" parameter

    all_pts = verts.union(set(T.pts))
    min_x = min(p.x() for p in all_pts)
    max_x = max(p.x() for p in all_pts)
    min_y = min(p.y() for p in all_pts)
    max_y = max(p.y() for p in all_pts)

    dx = (max_x-min_x)*(margin-1)/2
    dy = (max_y-min_y)*(margin-1)/2

    sw = Point(min_x-dx,min_y-dy)
    nw = Point(min_x-dx,max_y+dy)
    se = Point(max_x+dx,min_y-dy)
    ne = Point(max_x+dx,max_y+dy)

    walls = [
        Segment(se,ne),
        Segment(ne,nw),
        Segment(nw,sw),
        Segment(sw,se),
    ]

    for tri in tris.keys():

        for a,b,c in tri.adj():
            if Segment(a,b) in T.hull_edges:

                if cw(a,b,c):
                    a,b = b,a
                
                # ccw(a,b,c): c lies left of a -> b
                sab = Segment(a,b)
                center = tri.circum().get_center()
                bisect = Segment(a,b).bisector()
                midp = sab.midpoint()
                endpt = None

                # hack to represent semi-infinite edges as long segments
                #   since matplotlib does not support semi-infinite lines
                for wall in walls:
                    inter = wall.intersect_line(bisect)

                    if inter is None:
                        continue
                    elif ccw(a,b,inter):
                        continue
                    else:
                        if endpt is None:
                            endpt = inter
                        elif distance_to(midp, inter) < distance_to(midp, endpt):
                            endpt = inter

                assert(endpt is not None)
                semi_infinite_edges.add(Segment(center, endpt))

    return verts, bounded_edges, semi_infinite_edges, walls

def draw_voronoi(pts, margin=2):
    T = Triangulation(pts, True, True)
    T.random_incremental()
    assert(T.validate())
    T.draw()

    verts, bounded_edges, semi_infinite_edges, bounding_box = get_voronoi(T, margin)
    
    for e in bounded_edges:
        e.draw(color='firebrick')

    for e in semi_infinite_edges:
        e.draw(color='firebrick', arrow=True)
    
    for v in verts:
        v.draw(color='firebrick')

    for p in T.pts:
        p.draw(color='black')
    
    for w in bounding_box:
        w.draw(color='black')

    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    random.seed(300)
    P = sample_integer_points(10)
    draw_voronoi(P)