"""
CTSM-specific test that first runs the subset_data region tool and then ensures
that CTSM does not fail using the just-generated input files
"""

from subsetdata import SUBSETDATASHARED


class SUBSETDATAREGION(SUBSETDATASHARED):
    def __init__(self, case):
        """
        initialize an object interface to the SMS system test
        """

        lat1 = -9
        lat2 = -7
        lon1 = 291
        lon2 = 293

        # Don't need to include things that are added during SUBSETDATASHARED.__init__()
        subset_data_cmd = [
            "tools/site_and_regional/subset_data",
            "region",
            "--lat1",
            str(lat1),
            "--lat2",
            str(lat2),
            "--lon1",
            str(lon1),
            "--lon2",
            str(lon2),
            "--create-mesh",
            "--create-domain",
            "--create-surface",
            "--crop",
            "--create-landuse",
            "--surf-year",
            "1850",
        ]

        super().__init__(case, subset_data_cmd)
