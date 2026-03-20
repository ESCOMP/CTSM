"""
Constants used in crop calendar scripts
"""

# Define conversion multipliers, {from: {to1, to2, ...}, ...}
multiplier_dict = {
    # Mass
    "g": {
        "Mt": 1e-12,
    },
    "t": {
        "Mt": 1e-6,
    },
    # Volume
    "m3": {
        "km3": 1e-9,
    },
    # Yield
    "g/m2": {
        "t/ha": 1e-6 * 1e4,
    },
}

# Minimum harvest threshold allowed in PlantCrop()
# Was 50 before cropcal runs 2023-01-28
DEFAULT_GDD_MIN = 1.0
