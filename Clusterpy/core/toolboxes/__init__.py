import os 
import sys
from Clusterpy.core.toolboxes.cluster import *

folders = os.listdir(os.path.split(__file__)[0])
if "rimaps" in folders:
    try:
        from rimaps import *
        rimapsActive = True
    except Exception as e:
        print("Some functions are not available, reason:", e)