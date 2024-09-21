from Clusterpy.core import * 

def CPhelp(function = ''):
    """ClusterPy's official help system

    **Callable functions are:**

        * new: creates an empty Layer.
        * load: load a ClusterPy project (<file>.CP).
        * importArcData: creates a new Layer from a shapefile (<file>.shp).
        * createPoints: creates a new Layer with points uniformly distributed on space.
        * createGrid: creates a new Layer with a regular lattice.
        * importShape: reads the geographic information stored on a shapefile.

    For more information about any function just type ''CPhelp('<function>')'' 
    or read the official documentation available on 'documentation <www.rise-group.org>'
    
    **Examples**

    To see the help of a class, in this case ''Layer'', type:
    
    >>> import clusterpy
    >>> clusterpy.CPhelp("Layer")

    For a specific function, just type the name of the function:

    **Example 1**
    
    >>> import clusterpy
    >>> clusterpy.CPhelp("importArcData")

    **Example 2**

    >>> import clusterpy
    >>> clusterpy.CPhelp("new")
    """
    if not function:
        print(CPhelp.__doc__)
    else:
        try:
            exec ('print('+ function + '.__doc__)')
        except:
            print ("Invalid Function, to see available functions execute \
            'CPhelp()'")