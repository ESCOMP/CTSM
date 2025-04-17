"""General-purpose utility functions"""

import logging

import requests

from ctsm.utils import abort

logger = logging.getLogger(__name__)


def download_file(url, fname, timeout=30):
    """
    Function to download a file.
    Args:
        url (str):
            url of the file for downloading
        fname (str) :
            file name to save the downloaded file.
        timeout (number, optional) :
            time in seconds to wait for response before exiting

    Raises:
        Error :
            When the file is not available on the server (status_code:404)
        Error:
            When download fails for any reason.
    """
    try:
        response = requests.get(url, timeout=timeout)

    # pylint: disable=broad-except
    except Exception as err:
        logger.warning("The server could not fulfill the request.")
        logger.warning("Something went wrong in downloading: %s", fname)
        err_msg = "Couldn't download file " + fname + "-- Error code:" + err
        abort(err_msg)

    with open(fname, "wb") as this_f:
        this_f.write(response.content)

    # -- Check if download status_code
    if response.status_code == 200:
        logger.info("Download finished successfully for : %s", fname)

    elif response.status_code == 404:
        logger.warning("This file is not available on the server: %s", fname)
        err_msg = "Couldn't download file " + fname + "-- Error code: " + "404"
        abort(err_msg)
