"""
CTSM-specific test that first runs the subset_data point tool and then ensures
that CTSM does not fail using the just-generated input files
"""

from subsetdata import SUBSETDATASHARED


class SUBSETDATAPOINT(SUBSETDATASHARED):
    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """

        lat = 45.402252
        lon = -92.798085

        # Don't need to include things that are added during SUBSETDATASHARED.__init__()
        subset_data_cmd = [
            "tools/site_and_regional/subset_data",
            "point",
            "--lat",
            str(lat),
            "--lon",
            str(lon),
            "--create-surface",
            "--crop",
            "--create-landuse",
            "--surf-year",
            "1850",
            "--create-datm",
            "--datm-syr",
            "1901",
            "--datm-eyr",
            "1901",
        ]

        super().__init__(case, subset_data_cmd)
