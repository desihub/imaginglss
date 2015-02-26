# A simple pure-Python class to implement a KD tree
# and basic search for lists of n-dimensional data.
#
# Algorithm based on article in Wikipedia.
# http://en.wikipedia.org/wiki/Kd-tree



class KDtree:
    """
    KDtree:
    A pure-Python class to do some simple KD tree construction and
    search.  Everything is manipulated through lists.
    """
    __author__ = "Martin White"
    __version__ = "1.0"
    __email__  = "mwhite@berkeley.edu"
    class TreeNode:
        """
        An empty class for the nodes.
        """
        pass
    def kdtree(self,point_list,split_dim=0):
        """
        kdtree(self,point_list,split_dim=0):
        Does the work of constructing the KD tree in terms of "TreeNode"s.
        """
        if not point_list:
            return(None)
        axis = split_dim%(self.Ndim)
        # Sort point list and choose median as pivot element
        point_list.sort(key=lambda point: point[axis])
        median = len(point_list) // 2 # choose median
        # Create node and construct subtrees
        node = self.TreeNode()
        node.location   = point_list[median]
        node.left_child = \
          self.kdtree(point_list[:median],(split_dim+1)%self.Ndim)
        node.right_child= \
          self.kdtree(point_list[median+1:],(split_dim+1)%self.Ndim)
        return(node)
    def treelist(self,node):
        """
        treelist(self,node):
        Returns a list of all points in the tree (or subtree).
        """
        if node == None:
            raise ValueError,"node is None in treelist"
        pl = [ node.location ]
        if node.left_child != None:
            pl +=  self.treelist(node.left_child)
        if node.right_child != None:
            pl += self.treelist(node.right_child)
        return(pl)
    def __init__(self, point_list, split_dim=0):
        """
        __init__(self, point_list, split_dim=0):
        Instantiates the class, using kdtree above.
        """
        self.Ndim = len(point_list[0]) # Assume all the same dimension
        self.node = self.kdtree(point_list)
    def simplesearch(self,point_value,depth):
        """
        simplesearch(self,point_value,depth):
        Walks down the tree "depth" steps, then returns all points in
        that subtree.  If depth is a mutiple of the number of dimensions
        then you get a spatially compact region.
        """
        # Walk down "depth" levels.
        idim = 0
        cur  = self.node
        for i in range(depth):
            if point_value[idim]<=cur.location[idim]:
                if cur.left_child != None:
                    cur = cur.left_child
                else:
                    break
            else:
                if cur.right_child != None:
                    cur = cur.right_child
                else:
                    break
            idim = (idim+1)%self.Ndim
        # Return a list of all points in the tree "below" this.
        return(self.treelist(cur))
