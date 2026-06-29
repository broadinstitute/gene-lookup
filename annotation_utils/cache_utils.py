import functools
import gzip
import hashlib
import json
import os
import pandas as pd
import re
import time

CACHE_DIR = os.path.expanduser("~/.annotations")

# Read FORCE_DOWNLOAD at call time (not import time) so callers can set the env var after
# importing modules that wrap themselves with these decorators.
def _force_download():
    return os.getenv("FORCE_DOWNLOAD") == "1"


def _cache_file_path(function_name, args, kwargs, ext="tsv.gz"):
    """Return the cache file path for a decorated function called with the given args/kwargs."""
    h = hashlib.sha256(f"{function_name} {args} {frozenset(sorted(kwargs.items()))}".encode()).hexdigest()[:10]
    return os.path.join(CACHE_DIR, re.sub("^get_", "", function_name) + f".{h}.{ext}")


def read_cached_table(function_name, *args, **kwargs):
    """Return the cached DataFrame for the given cache_data_table function regardless of age.

    Intended as a fallback when a fresh download fails: it returns the last successfully cached
    table even if it's older than the normal 5-day freshness window, or None if no cache exists.
    """
    cache_file_path = _cache_file_path(function_name, args, kwargs)
    if os.path.isfile(cache_file_path):
        return pd.read_table(cache_file_path)
    return None


def cache_data_table(get_table_func):
    """Decorator that caches the pandas DataFrame returned by the decorated function.
    It's intended for functions that take a relatively long time to retrieve some table over the network.
    Before calling the decorated function, the decorator checks whether result already exists in the
    cache dir (~/.annotations). If yes, it just reads the table from disk and returns it.
    If no, it calls the function and then saves the result table to ~/.annotations before returning it.
    """

    @functools.wraps(get_table_func)
    def wrapper(*args, **kwargs):
        # create cache dir
        if not os.path.isdir(CACHE_DIR):
            os.mkdir(CACHE_DIR)

        # check if cached file already exists
        cache_file_path = _cache_file_path(get_table_func.__name__, args, kwargs)

        # use the cached file if it's less than 5 days old
        if (
            not _force_download() and
            os.path.isfile(cache_file_path) and
            os.path.getmtime(cache_file_path) > time.time() - 5 * 24 * 60 * 60
        ):
            df = pd.read_table(cache_file_path)
            print(f"Read {len(df):,d} rows from cache file {cache_file_path}")
            return df

        # call the underlying function
        df = get_table_func(*args, **kwargs)

        # save result to cache
        df.to_csv(cache_file_path, header=True, index=False, sep="\t")

        print(f"Saved {len(df):,d} rows to cache file {cache_file_path}")
        return df

    return wrapper


def cache_json(get_json_func):
    """Decorator that caches the json returned by the decorated function.
    It's intended for functions that take a relatively long time to retrieve some json over the network.
    Before calling the decorated function, the decorator checks whether result already exists in the
    cache dir (~/.annotations). If yes, it just reads the json from disk and returns it.
    If no, it calls the function and then saves the result json to ~/.annotations before returning it.
    """

    @functools.wraps(get_json_func)
    def wrapper(*args, **kwargs):

        # create cache dir
        if not os.path.isdir(CACHE_DIR):
            os.mkdir(CACHE_DIR)

        # check if cached file already exists
        cache_file_path = _cache_file_path(get_json_func.__name__, args, kwargs, ext="json.gz")

        # use the cached file if it's less than 5 days old
        if (
            not _force_download() and
            os.path.isfile(cache_file_path) and
            os.path.getmtime(cache_file_path) > time.time() - 5 * 24 * 60 * 60
        ):
            return json.load(gzip.open(cache_file_path, "rt"))

        # call the underlying function
        json_data = get_json_func(*args, **kwargs)

        # save result to cache
        with gzip.open(cache_file_path, "wt") as f:
            json.dump(json_data, f, indent=2)

        return json_data

    return wrapper