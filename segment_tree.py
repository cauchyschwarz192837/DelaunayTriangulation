from primitives import *

class SegmentTreeAuxSet(): # slower because just using set and not balanced AVL
    '''a class to be used at every node of a SegmentTree to store its segments.

        (In lecture we presented a more involved one backed by an AVL tree, whereas
        this one uses a simple Python set and thus will have worse asymptotic performance.)'''
    
    def __init__(self, interval):
        '''intitialize an empty segment set for an interval tree
            whose node has the given interval'''
        self.segs = set()
        self.interval = interval

    def insert(self, seg):
        '''inserts this segment to this set'''
        self.segs.add(seg)

    def delete(self, seg):
        '''deletes the segment from this set'''
        self.segs.remove(seg)

    def get_segs(self):
        '''return a set of this set's segments'''
        return self.segs

    def vertical_shoot(self, p : Point):
        '''given a point p in self.interval, return the pair (above_seg, above_point) where
        `above_seg` is the lowest Segment in self.segs visible upwards from p (possibly containing p),
        and above_point is the point on `above_seg` visible from p. If there is no such segment above or
        containing p, return (None, None).
        
        Target runtime is O(1) per Segment in self.segs. You do NOT need to implement the faster logarithmic method.'''

        assert(self.interval.contains_1d_point(p.x_proj())) # verify the given point's x-coordinate lies in this set's assigned interval

        above_seg = None
        above_point = None
        
        vertical_line = p.vertical_line_thru()

        for s in self.segs:
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

        # TODO: Implement this method (Task 3)

        return (above_seg, above_point)

class SegmentTree():

    @classmethod
    def from_2d_points(cls, points):
        '''return a SegmentTree built over the x-coordinates of the given points'''

        x_coords = sorted(set(p.x_proj() for p in points))
        return cls(x_coords)

    def __init__(self, x_coords):
        '''build an empty segment tree on the given list 1d points, x_coords'''
        self.pts = x_coords
        self.interval = Interval(x_coords[0], x_coords[-1]) # big guy
        self.aux = SegmentTreeAuxSet(self.interval)
        
        if len(x_coords) == 2: # base case
            self.left = None
            self.right = None
            self.split = None
            return
        
        n = len(x_coords)
        l_coords = x_coords[:n//2+1] # median point belongs to both sides
        r_coords = x_coords[n//2:]
        self.split = l_coords[-1]
        
        self.left = SegmentTree(l_coords)
        self.right = SegmentTree(r_coords)


    def insert(self, seg : Segment):
        '''Given a non-vertical segment, insert it into the auxiliary data structure at
        its canonical nodes in the SegmentTree.'''

        assert not seg.is_vertical(), "Only non-vertical segments can be inserted."

        seg_extent = seg.x_extent()
        if seg_extent.contains(self.interval):
            self.aux.insert(seg)
            return

        if self.left and seg_extent.intersects(self.left.interval):
            self.left.insert(seg)
        if self.right and seg_extent.intersects(self.right.interval):
            self.right.insert(seg)

        # TODO: Complete for Task 3


    def delete(self, seg : Segment):
        '''given a non-vertical segment, delete it from the auxiliary data structure at
        its canonical nodes in this SegmentTree.
        
        Hint: "Undo" the insert() method.'''

        seg_extent = seg.x_extent()
        if seg_extent.contains(self.interval):
            self.aux.delete(seg)
            return
        
        if self.left and seg_extent.intersects(self.left.interval):
            self.left.delete(seg)
        if self.right and seg_extent.intersects(self.right.interval):
            self.right.delete(seg)

        # TODO: Complete for Task 3

    def vertical_shoot(self, p : Point):
        '''given a point whose x-coordinate lies in this node's interval,
        return the segment of this triangulation that either
        contains this point OR the lowest segment visible upwards from p, and the
        visible point on that segment. If the segment contains p, the visible point
        is p itself.
        
        Implementation hints:
            - For a given Segment, its left and right (2D Point) endpoints can be accessed with .left and .right.
            - The 1D projection of a 2D Point p (i.e., a 1D point with its x-coordinate) can be obtained with p.x_proj(),
                which is a OneDPoint object. See these classes and their methods in primitives.py.
            - Two OneDPoint objects' x-coordinates can be compared using <, <=, ==, >, >=.'''

        above_seg, above_point = None, None

        retminim = float('inf')
        vertical_line = p.vertical_line_thru()
        for s in self.aux.get_segs():
            if s.contains_point(p):
                above_seg = s
                above_point = p
            intersection = vertical_line.intersect_segment(s)
            if intersection:
                if intersection.y() > p.y() and intersection.y() < retminim:
                    retminim = intersection.y()
                    above_seg = s
                    above_point = intersection

        if self.left and p.x_proj() <= self.split:
            left_seg, left_point = self.left.vertical_shoot(p)
            if left_seg and left_point.y() < retminim:
                above_seg, above_point = left_seg, left_point
                retminim = left_point.y()

        if self.right and p.x_proj() >= self.split:
            right_seg, right_point = self.right.vertical_shoot(p)
            if right_seg and right_point.y() < retminim:
                above_seg, above_point = right_seg, right_point

        # TODO: Complete for Task 3

        return (above_seg, above_point)
    
    def stabbing_query(self, q):
        '''returns the set of all Segments in this SegmentTree `stabbed` by the
        vertical line through the given query point q; i.e., all Segments
        that contain q's x-coordinate.'''

        ret = set()

        if q.x_proj() < self.interval.left or q.x_proj() > self.interval.right:
            return ret
        
        ret = self.aux.get_segs()

        if self.left:
            ret = ret.union(self.left.stabbing_query(q))

        if self.right:
            ret = ret.union(self.right.stabbing_query(q))
        
        return ret

    def gather(self):
        '''return a set of all segments stored in this subtree'''
        ret = self.aux.get_segs()
        
        if self.left:
            ret = ret.union(self.left.gather())
        
        if self.right:
            ret = ret.union(self.right.gather())

        return ret
        
    def draw_stabbing_query(self, q):
        all_segs = self.gather()
        stabbed = self.stabbing_query(q)

        for s in all_segs:
            s.draw()
        for s in all_segs.difference(stabbed):
            s.draw(color='gray')

        for s in stabbed:
            s.draw(color='red')

        q.vertical_line_thru().draw()

    def draw_vertical_shoot(self, q):
        above, above_point = self.vertical_shoot(q)
        all_segs = self.gather()
        
        for s in all_segs:
            s.draw()

        above.draw(color='forestgreen')
        q.vertical_line_thru().draw()
        above_point.draw(color='red')
        q.draw(color='purple')

    def draw(self, depth=0, top=None):
        COLORS = ['#999999', '#F781BE', '#E31929', '#3A7FB5', '#4CAF52', '#984FA0', '#FE7E2A', '#A55630']

        if top is None:
            all_segs = self.gather()
            retminim = min(min(s.p1.y(), s.p2.y()) for s in all_segs)
            top = retminim-10
            for p in self.pts:
                p.lift(0).vertical_line_thru().draw(color='lightgray',dashed=True)

        # get vertical lines through interval endpoints
        left_vert = self.interval.left.lift(0).vertical_line_thru()
        right_vert = self.interval.right.lift(0).vertical_line_thru()
        if self.split:
            split_pt = self.split.lift(top-10*depth)
            Segment(split_pt.translate(0,-10), split_pt).draw(color='gray',arrow=True)

        # for each segment stored at this node, color the portion whose extent lies in this interval
        for s in self.aux.get_segs():
            clip_left = s.intersect_line(left_vert)
            clip_right = s.intersect_line(right_vert)
            Segment(clip_left, clip_right).draw(color=COLORS[depth % len(COLORS)])
            clip_left.draw(color='gray')
            clip_right.draw(color='gray')
        
        # draw this interval at y-coordinate shifted down based on depth
        seg = Segment(self.interval.left.lift(top-10*depth), self.interval.right.lift(top-10*depth))
        seg.draw(color=COLORS[depth % len(COLORS)])

        if self.left:
            self.left.draw(depth+1, top)            

        if self.right:
            self.right.draw(depth+1, top)

if __name__=='__main__':
    from delaunay import *
    import random

    random.seed(290)

    P = sample_integer_points(15)
    T = Triangulation(P, False, False)
    # T.random_incremental()
    # T.naive_delaunay()
    T.draw()
    tree = SegmentTree.from_2d_points(P)
    for seg in T.edges:
        tree.insert(seg)
    
    tree.draw()
    plt.show()

    tree.draw_stabbing_query(Point(35,15))
    plt.show()

    tree.draw_vertical_shoot(Point(35,15))
    plt.show()
