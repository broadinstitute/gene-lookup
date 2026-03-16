import functools
import gzip
import hashlib
import json
import os
import pandas as pd
import re
import time

CACHE_DIR = os.path.expanduser("~/.annotations")

FORCE_DOWNLOAD = os.getenv("FORCE_DOWNLOAD") == "1"

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
        function_name = get_table_func.__name__
        h = hashlib.sha256(f"{function_name} {args} {frozenset(sorted(kwargs.items()))}".encode()).hexdigest()
        h = h[:10]
        filename = re.sub("^get_", "", function_name) + f".{h}.tsv.gz"
        cache_file_path = os.path.join(CACHE_DIR, filename)

        # use the cached file if it's less than 5 days old
        if (
            not FORCE_DOWNLOAD and
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
        function_name = get_json_func.__name__
        h = hashlib.sha256(f"{function_name} {args} {frozenset(sorted(kwargs.items()))}".encode()).hexdigest()
        h = h[:10]
        filename = re.sub("^get_", "", function_name) + f".{h}.json.gz"
        cache_file_path = os.path.join(CACHE_DIR, filename)

        # use the cached file if it's less than 5 days old
        if (
            not FORCE_DOWNLOAD and
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