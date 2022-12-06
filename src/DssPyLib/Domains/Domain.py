class Domain:
    """
    A class to represent a geometry.
    
    ...

    Attributes
    ----------
    width : float
        width of the geometry
    height : float
        height of the geometry
    position : list
        position of the most significant point of the geometry
    relativeMeshDensity : float
        relative mesh density with respect to the global geometry
    angle : float
        angle of rotation of the geometry around the most significant point of the geometry
    isChild : boolean
        states if a geometry is a child of another geometry
    isParent : boolean
        states if a geometry is a parent of another geometry
    childDomains : list
        contains list of child domains in the geometry
    parentDomains : int
        contains the number of parent domain of current geometry
    non : int
        number of nodes to be forced manually along the width of the most significant perimeter

    Methods
    -------
    calculateRelativeWidth(gWidth):
        calculates relative width of a geometry with respects to the global geometry width.
    calculateWidthNodes(gNon):
        assigns number of nodes along the width of a geometry according to non or relative mesh density.
    """
    
    def __init__(self, width, height, position, relativeMeshDensity, angle, isChild, isParent, childDomains, parentDomains, non):
        """
        Constructs all the necessary attributes for the geometry object.

        Parameters
        ----------
        width : float
            width of the geometry
        height : float
            height of the geometry
        position : list
            position of the most significant point of the geometry
        relativeMeshDensity : float
            relative mesh density with respect to the global geometry
        angle : float
            angle of rotation of the geometry around the most significant point of the geometry
        isChild : boolean
            states if a geometry is a child of another geometry
        isParent : boolean
            states if a geometry is a parent of another geometry
        childDomains : list
            contains list of child domains in the geometry
        parentDomains : int
            contains the number of parent domain of current geometry
        non : int
            number of nodes to be forced manually along the width of the most significant perimeter
        """

        self.width = width
        self.height = height
        self.position = position
        self.relativeMeshDensity = relativeMeshDensity
        self.angle = angle
        self.isChild = isChild
        self.isParent = isParent
        self.childDomains = childDomains
        self.parentDomains = parentDomains
        self.non = non

    def calculateRelativeWidth(self, gWidth):
        """
        calculates relative width of a geometry with respects to the global geometry width.

        Parameters
        ----------
        gwidth : float
            width og gloabl geometry

        Returns
        -------
        None
        """

        self.relativeWidth = self.width/gWidth

    def calculateWidthNodes(self, gNon):
        """
        assigns number of nodes along the width of a geometry according to non or relative mesh density.

        Parameters
        ----------
        gNon : int
            number of nodes on the width of global boundary

        Returns
        -------
        None
        """
        
        if self.non == 0:
            self.widthNodes = round(self.relativeWidth*self.relativeMeshDensity[0]*gNon+1)
        else:
            self.widthNodes = self.non

    
