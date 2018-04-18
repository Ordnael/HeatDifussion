# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 23:58:05 2018

Library script containing the materials to analyze.

@author: Leandro
"""

"""
Import statements
"""

import numpy as np


"""
Dictionary of materials.
"""

# Water.
water = {}
water['temp'] = np.array([1.0, 1.5, 2.0, 2.5])
water['coef'] = np.array([1.0, 1.1375, 1.35, 1.6375])
water['degr'] = 2

