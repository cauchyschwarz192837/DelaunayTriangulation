from primitives import *
import random
import matplotlib.pyplot as plt
from itertools import islice
from graham import graham
import numpy as np

from segment_tree import *
from segment_tree import SegmentTree

class Triangulation():

    def __init__(self, pts, use_tree=False, make_legal=False, DRAW=False, SAVE_TO_GIF=False):
        '''Create a new Triangulation given a list of 2D Points by first inserting all points
        on its convex hull (whose edges must be in the triangulation), then inserts the rest in some sorted order.
        
        Optional parameters:
        - If `use_tree` is True, then a SegmentTree will be used to identify the segment above a given point to insert.
        
        - If `make_legal` is True, then the legalize() method is called as new segments of the triangulation are created.
        
        - If `SAVE_TO_GIF` is True, then the .save_plot() method will save the figure currently in `plt` as one frame of the GIF to be exported.
            Call self.show_plot() to either display the plot (if SAVE_TO_GIF is False) or save the current plot as one frame of the GIF to be exported.

        - If `DRAW` is True, self.DRAW can be used in the bodies of the methods to optionally draw for debugging/animation purposes.
        
        Attributes:
            pts
            hull_pts
            hull_edges
            edges
            adj
            tree
            make_legal
        '''
        self.DRAW = DRAW
        self.SAVE_TO_GIF = SAVE_TO_GIF
        self.TITLE = None

        self.pts = list(pts)
        self.verts = []
        
        # compute convex hull of pts
        hull = graham(pts)
        h = len(hull)
        self.hull_pts = set(hull)
        self.hull_edges = set(Segment(hull[i],hull[(i+1)%h]) for i in range(h))

        # initialize empty edges list and adjacency map
        self.edges = []
        self.adj = {}

        if use_tree:
            self.tree = SegmentTree.from_2d_points(pts)
        else:
            self.tree = None

        self.make_legal = make_legal

        # initialize this triangulation as a triangle with first 3 points on the hull
        p1,p2,p3 = hull[:3]
        self.add_segment(p1, p2)
        self.add_segment(p2, p3)
        self.add_segment(p3, p1)

        if self.DRAW:
            self.draw()
            self.show_plot()

        # for the remaining points on the hull, add them one-by-one by connecting to
        #   the first and previous points on the hull
        for i in range(3,h):
            self.add_segment(hull[i], hull[0]) # new outer edge
            self.add_segment(hull[i], hull[i-1]) # new inner edge

            if self.DRAW:
                self.draw()
                Segment(hull[i], hull[0]).draw(color='black')
                Segment(hull[i], hull[i-1]).draw(color='black')
                self.show_plot()
            
            if self.make_legal: 
                self.legalize(hull[i], hull[i-1], hull[0]) # legalize inner edge
    
    def random_incremental(self):
        '''add the rest of the points in self.pts to the tree (assumes the hull points were
        added during the constructor, __init__()).'''
        
        random.shuffle(self.pts)    
        for p in self.pts:
            if p not in self.hull_pts: # skip points already accounted for
                self.insert_point(p)

    def naive_ray_shoot(self, p : Point):
        '''given a point p, return the segment of this triangulation that either
        contains this point OR the lowest segment visible upwards from p, and the
        visible point on that segment. If the segment contains p, the visible point
        is p itself.

        ASSUMPTION: no vertex of the current triangulation has the same x-coordinate as p'''        

        above_seg = None
        above_point = None

        vertical_line = p.vertical_line_thru()

        for s in self.edges:
            if s.contains_point(p):
                above_seg = s
                above_point = p
                break
            else:
                inters = vertical_line.intersect_segment(s)
                if inters is not None:
                    if s.contains_point(inters):
                        if inters.is_above(p) or inters.equal_y(p):
                            if above_point is None or inters.is_below(above_point):
                                above_seg = s
                                above_point = inters

        # TODO: Implement this method (Task 1)
        # Hint: A vertical Line object through a given (2D) Point object "p"
        #           can be obtained with p.vertical_line_thru()

        return (above_seg, above_point)

    def naive_delaunay(self):
        '''While there are illegal edges in the triangulation, flip them.
        When it terminates, this triangulation is Delaunay.

        ASSUMPTION: all points of self.pts have been inserted into the triangulation'''

        stack = set()  # Use a set to avoid duplicates
        for s in self.edges:
            stack.add((s.p1, s.p2))

        while len(stack) > 0:
            a, b = stack.pop()
            res = self.is_illegal(a, b)

            if res is (None, None):
                continue
            else:
                c, d = res
                if Segment(a, b) in self.edges:
                    self.remove_segment(a, b)
                if Segment(b, a) in self.edges:
                    self.remove_segment(b, a)
                if Segment(c, d) not in self.edges and Segment(d, c) not in self.edges:
                    self.add_segment(c, d) # do not need to add to stack, must be legal

                '''add the given segment to this triangulation, updating its adjacency map and segment tree
                seg = Segment(a,b)
                self.edges.append(seg)'''

                edges_to_check = [
                    (a, c), (c, b), (d, b), (d, a)
                ]
                for edge in edges_to_check:
                    if Segment(edge[0], edge[1]) in self.edges or Segment(edge[1], edge[0]) in self.edges:
                        stack.add(edge)
            
            # TODO: Finish this implementation (Task 1):
            #  - Call .is_illegal() to check the popped edge ab
            #  - If necessary, flip this edge by calling self.remove_segment() and self.add_segment() appropriately
            #  - Finally, add appropriate segments onto the "stack"



    def legalize(self, p, a, b):
        '''Given a newly-inserted point p and edge (a, b), this method checks if edge (a, b) is illegal.
        If it is, it flips the edge to edge (p, q), where q is the point opposite to edge (a, b) in the
        adjacent triangle. It then recursively legalizes affected edges.'''

        seg = Segment(a, b)
        assert(seg in self.edges)
        if seg in self.hull_edges:
            return

        incident_a = self.get_incident(a) # CW order
        idx_b = incident_a.index(b)

        idx_temp = (idx_b + 1) % len(incident_a) # need to mod
        temp = incident_a[idx_temp]
        while temp not in self.adj[b]: # this loop is not necessary but can't prove it yet so leave it first
            idx_temp = (idx_temp + 1) % len(incident_a)
            temp = incident_a[idx_temp]
        idx_q = idx_temp
        q = temp

        if ccw(a, p, b):
            ccw_points = [p, b, q, a]
        else:
            ccw_points = [p, a, q, b]

        if q == p:
            idx_temp = (idx_b - 1) % len(incident_a) # need to mod
            temp = incident_a[idx_temp]
            while temp not in self.adj[b]:
                idx_temp = (idx_temp - 1) % len(incident_a)
                temp = incident_a[idx_temp]
            idx_q = idx_temp
            q = temp

        if not self.verify_convex(ccw_points):
            return

        circum = Circle(p, a, b)
        if circum.in_circle(q):
            self.remove_segment(a, b)
            self.add_segment(p, q)

            self.legalize(p, q, a)
            self.legalize(p, q, b)
        else:
            return

        # TODO: Complete this method (Task 2)
        #   1) Find the vertex of the triangle using segment ab which is not p
        #   2) Check if the quadrilateral of the four points are convex
        #   3) If so, determine if ab is illegal, and if yes, flip it and finish accordingly.

        # You may find the primitives.Circle.in_circle() and Triangulation.verify_convex() methods useful.
    

    def is_illegal(self, a, b):
        '''return (None, None) if segment ab is legal, otherwise return the two points c,d on the 
        convex quadrilateral with a,b for which segment ab is illegal and cd is legal.'''
        
        if not Segment(a,b) in self.edges or Segment(a,b) in self.hull_edges:
            return (None, None)

        c = self.get_ccw_neighbor(a,b)
        d = self.get_cw_neighbor(a,b)

        if not self.verify_convex([c,a,d,b]):
            return (None, None)
        
        if self.DRAW:
            Triangle(a,b,c).draw(color='lightgray')
            Triangle(a,b,d).draw(color='lightgray')
            self.draw()
            Segment(a,b).draw(color='red')
            a.draw(color='darkorange')
            b.draw(color='darkorange')
            c.draw(color='darkorange')
            self.show_plot()

        circ = Circle(a,b,c)

        if self.DRAW:
            Triangle(a,b,c).draw(color='lightgray')
            Triangle(a,b,d).draw(color='lightgray')
            self.draw()
            circ.draw()
            Segment(a,b).draw(color='red')
            a.draw(color='darkorange')
            b.draw(color='darkorange')
            c.draw(color='darkorange')
            self.show_plot()
            
        if not circ.in_circle(d):
            return (None, None)
        
        if self.DRAW:
            Triangle(a,b,c).draw(color='lightgray')
            Triangle(a,b,d).draw(color='lightgray')
            self.draw()
            circ.draw()
            Segment(a,b).draw(color='red')
            Segment(c,d).draw(color='green')
            a.draw(color='darkorange')
            b.draw(color='darkorange')
            c.draw(color='darkorange')
            self.show_plot()
            
        return (c,d)

    def verify_convex(self, pts):
        '''given list of points in CCW order, return True if and only if polygon is convex'''
        return all(not cw(pts[(i-1)%len(pts)],pts[i],pts[(i+1)%len(pts)]) for i in range(len(pts)))

    def add_segment(self, a, b):
        '''add the given segment to this triangulation, updating its adjacency map and segment tree'''

        seg = Segment(a,b)

        self.edges.append(seg)
        
        if seg.p1 not in self.adj:
            self.adj[seg.p1] = set()
            
        if seg.p2 not in self.adj:
            self.adj[seg.p2] = set()

        self.adj[seg.p1].add(seg.p2)
        self.adj[seg.p2].add(seg.p1)

        if self.tree and not seg.is_vertical():
            self.tree.insert(seg)

    def remove_segment(self, a, b):
        '''remove the given segment from this triangulation, updating its adjacency map and segment tree'''
        
        seg = Segment(a,b)

        self.edges.remove(seg)

        self.adj[seg.p1].remove(seg.p2)
        self.adj[seg.p2].remove(seg.p1)

        if self.tree and not seg.is_vertical():
            self.tree.delete(seg)

    def get_incident(self, p):
        '''given a point p of the triangulation, return a sorted list of its of adjacent points
        in clockwise order'''

        adj = self.adj[p]
        adj = sorted( adj, key=lambda q: -p.angle(q) ) # sorts CW around p
        return adj
    
    def get_cw_neighbor(self, a, b):
        '''given adjacent points a,b of this triangulation,
        return the vertex adjacent to a in clockwise order after b'''

        adj = self.get_incident(a)
        idx = adj.index(b)
        return adj[(idx+1)%len(adj)]

    def get_ccw_neighbor(self, a, b):
        '''given adjacent points a,b of this triangulation,
        return the vertex adjacent to a in counter-clockwise order after b'''

        adj = self.get_incident(a)
        idx = adj.index(b)
        return adj[(idx-1)%len(adj)]
    
    def insert_point(self, p):
        '''given a point p, insert it to the triangulation then modify it into a valid
        triangulation. If use_tree=True use self.tree to find p's visible segment, otherwise
        use the naive method. If legalize=True, legalize all relevant segments recursively.

        ASSUMPTION: the given point p is contained in the interior of a triangle (face) of this triangulation
        OR it lies on the interior of a segment (edge) not on the convex hull'''

        if self.tree:
            above, above_point = self.tree.vertical_shoot(p)
        else:
            above, above_point = self.naive_ray_shoot(p)

        if self.DRAW:
            self.draw()
            above.draw(color='darkorange')
            above_point.draw(color='orange')
            if p != above_point:
                Segment(p,above_point).draw(color='black', arrow=True)
            p.draw(color='red')
            self.show_plot()

        a = above.left
        b = above.right

        c = self.get_ccw_neighbor(a,b)
        d = self.get_cw_neighbor(a,b)

        if above.contains_interior_point(p):

            self.remove_segment(a,b)

            if self.DRAW:
                self.draw()
                Segment(p,a).draw(color='darkorange')
                Segment(p,b).draw(color='darkorange')
                Segment(p,c).draw(color='darkorange')
                Segment(p,d).draw(color='darkorange')
                p.draw(color='red')
                self.show_plot()

            self.add_segment(p, a)
            self.add_segment(p, b)
            self.add_segment(p, c)
            self.add_segment(p, d)

            if self.make_legal:
                self.legalize(p, c, a)
                self.legalize(p, a, d)
                self.legalize(p, d, b)
                self.legalize(p, b, c)
            
            return

        # if p not on a segment, above or below is an edge    

        self.add_segment(p, a)
        self.add_segment(p, b)
        self.add_segment(p, d)
        
        if self.make_legal:
            self.legalize(p, b, a)
            self.legalize(p, a, d)
            self.legalize(p, d, b)

        if self.DRAW:
            self.draw()
            Segment(p,a).draw(color='firebrick')
            Segment(p,b).draw(color='firebrick')
            Segment(p,d).draw(color='firebrick')
            p.draw(color='red')
            self.show_plot()
    
    def draw(self):
        for s in self.edges:
            s.draw(color='gray')
        
        for p in self.pts:
        # for p in self.adj.keys():
            p.draw()
    
    def get_triangles(self):
        '''return the bounded triangles of this Triangulation'''
        tris = {}

        for seg in set(self.edges).difference(self.hull_edges):
            a, b = seg.p1, seg.p2

            c = self.get_cw_neighbor(b, a)
            d = self.get_ccw_neighbor(b, a)

            tri_c = Triangle(a,b,c)
            tri_d = Triangle(a,b,d)

            if tri_c not in tris:
                tris[tri_c] = []
            if tri_d not in tris:
                tris[tri_d] = []

            tris[tri_c].append(tri_d)
            tris[tri_d].append(tri_c)
        
        return tris

    def get_scipy_reference(self):
        import numpy as np
        from scipy.spatial import Delaunay

        pz = list(list(p.p()) for p in self.pts)
        pz = np.array(pz)
        tri = Delaunay(pz)
        edges = set()
        for a,b,c in tri.simplices:
            pa = Point(int(pz[a][0]), int(pz[a][1]))
            pb = Point(int(pz[b][0]), int(pz[b][1]))
            pc = Point(int(pz[c][0]), int(pz[c][1]))
            sa = Segment(pa, pb)
            sb = Segment(pc, pb)
            sc = Segment(pa, pc)
            edges.add(sa)
            edges.add(sb)
            edges.add(sc)

        return edges
    
    def validate(self):
        '''using scipy, use their implementation then compare its triangulation's segments to ours'''
        scipy_edges = self.get_scipy_reference()
        return set(self.edges) == scipy_edges
    
    def show_plot(self,margin=1.05):
        '''outputs the current figure in plt to as a .png in the list self.FRAMES, otherwise displays as normal'''
        try:
            self.min_x
        except:
            self.FRAMES = []
            self.min_x = min(p.x() for p in self.hull_pts)
            self.max_x = max(p.x() for p in self.hull_pts)
            self.min_y = min(p.y() for p in self.hull_pts)
            self.max_y = max(p.y() for p in self.hull_pts)
            self.dx = (margin-1)*(self.max_x - self.min_x)
            self.dy = (margin-1)*(self.max_y - self.min_y)

        plt.axis([self.min_x-self.dx, self.max_x+self.dx, self.min_y-self.dy, self.max_y+self.dy])
        plt.gca().set_aspect('equal')
        if self.SAVE_TO_GIF:
            from PIL import Image # pip install pillow
            import io

            buf = io.BytesIO()
            plt.savefig(buf, format='png')
            plt.close()
            buf.seek(0)
            self.FRAMES.append(Image.open(buf))
        else:
            plt.show()

    def save_plot(self, fname, pause=1.0):
        '''saves the plots stored in self.FRAMES, populated by self.show_plot() with
        when SAVE_TO_GIF=True in the constructor. fname is the local file name, and pause
        is an optional parameter that is the number of seconds per frame (default 1.0).'''
        if self.SAVE_TO_GIF:
            self.FRAMES[0].save(
                fname,
                format='GIF',
                append_images=self.FRAMES[1:],
                save_all=True,
                duration=int(pause*1000),  # Time in milliseconds between frames
                loop=0         # 0 for infinite loop
            )

def sample_integer_points(n,xoffset=0,yoffset=0,sparsity=5):
    '''returns a set of points with distinct integer coordinates,
    as a result no two points lie on the same horizontal or vertical lines.'''
    xs = list(range(sparsity*n))
    ys = list(range(sparsity*n))
    random.shuffle(xs)
    random.shuffle(ys)
    return [Point(x+xoffset*sparsity,y+yoffset*sparsity) for x,y in islice(zip(xs,ys), n)]

if __name__=='__main__':
    random.seed(290)
    
    P = sample_integer_points(12)
    ###### Task 0: Getting Started ############################
    random.seed(290)
    T = Triangulation(P, use_tree=False, make_legal=False, DRAW=True, SAVE_TO_GIF=True) # sets up the Triangulation with its convex hull
    T.draw() # draw the final triangulation
    plt.show() # show it
    T.save_plot(fname='task0.gif')

    ###### Task 1: Create a Triangulation by naively flipping edges after all points have been added #######
    random.seed(290)
    T = Triangulation(P, use_tree=False, make_legal=False, DRAW=True, SAVE_TO_GIF=True) # sets up the Triangulation with its convex hull
    T.random_incremental() # adds the rest of the points to the Triangulation
    T.naive_delaunay() # while there is an illegal edge, flip it
    print(T.validate()) # print whether the triangulation is correct
    T.draw() # draw the final triangulation
    plt.show() # show it
    T.save_plot(fname='task1.gif')

    
    # ####### Task 2: Speed up Task 1 further by legalizing edges as they are created instead of at the end ######
    random.seed(290)
    T = Triangulation(P, use_tree=False, make_legal=True, DRAW=True, SAVE_TO_GIF=True) # sets up the Triangulation with its convex hull
    T.random_incremental() # adds the rest of the points to the Triangulation
    print(T.validate()) # print whether the triangulation is correct
    T.draw() # draw the final triangulation
    plt.show() # show it
    T.save_plot(fname='task2.gif', pause=0.8)
    
    
    ####### Task 3: Speed up Task 1 by using a SegmentTree for point location #######
    random.seed(290)
    T = Triangulation(P, use_tree=True, make_legal=True, DRAW=True, SAVE_TO_GIF=True) # sets up the Triangulation with its convex hull
    T.random_incremental() # adds the rest of the points to the Triangulation
    print(T.validate()) # print whether the triangulation is correct
    T.draw() # draw the final triangulation
    plt.show() # show it
    T.save_plot(fname='task3.gif')
    