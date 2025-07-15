"""
Constants and functions relating to PFTs
"""

MIN_PFT = 0  # bare ground
MIN_NAT_PFT = 1  # minimum natural pft (not including bare ground)
MAX_NAT_PFT = 14  # maximum natural pft
MAX_PFT_GENERICCROPS = 16  # for runs with generic crops
MAX_PFT_MANAGEDCROPS = 78  # for runs with explicit crops


def is_valid_pft(pft_num, managed_crops):
    """
    Given a number, check whether it represents a valid PFT (bare ground OK)
    """
    if managed_crops:
        max_allowed_pft = MAX_PFT_MANAGEDCROPS
    else:
        max_allowed_pft = MAX_PFT_GENERICCROPS

    return MIN_PFT <= pft_num <= max_allowed_pft
