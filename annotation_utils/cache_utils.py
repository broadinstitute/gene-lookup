import functools
import glob
import gzip
import hashlib
import json
import os
import pandas as pd
import re
import time

CACHE_DIR = os.path.expanduser("~/.annotations")
SOURCE_LAST_UPDATED_PATH = os.path.join(CACHE_DIR, "source_last_updated.json")

# Marks a DataFrame that was read from the cache rather than downloaded
FROM_CACHE_ATTR = "read_from_cache"

# Read FORCE_DOWNLOAD at call time (not import time) so callers can set the env var after
# importing modules that wrap themselves with these decorators.
def _force_download():
    return os.getenv("FORCE_DOWNLOAD") == "1"


def _cache_file_path(function_name, args, kwargs, ext="tsv.gz"):
    """Return the cache file path for a decorated function called with the given args/kwargs."""
    h = hashlib.sha256(f"{function_name} {args} {frozenset(sorted(kwargs.items()))}".encode()).hexdigest()[:10]
    return os.path.join(CACHE_DIR, re.sub("^get_", "", function_name) + f".{h}.{ext}")


def read_source_last_updated():
    """Return a {source_name: "YYYY-MM-DD"} dict recording when each source was last downloaded successfully."""
    if os.path.isfile(SOURCE_LAST_UPDATED_PATH):
        with open(SOURCE_LAST_UPDATED_PATH) as f:
            return json.load(f)
    return {}


def record_source_last_updated(source_name):
    """Record today's UTC date as the last successful download date for source_name.

    Sources whose download key expires (OMIM, dbNSFP) fall back to their stale cached table instead of
    failing the pipeline, so this is called only after a fresh download succeeds. The recorded date then
    stays put while a source is stuck on its stale cache, which makes that staleness visible on the website.
    """
    dates = read_source_last_updated()
    dates[source_name] = time.strftime("%Y-%m-%d", time.gmtime())
    with open(SOURCE_LAST_UPDATED_PATH, "wt") as f:
        json.dump(dates, f, indent=2, sort_keys=True)


def read_cached_table(function_name):
    """Return the newest cached DataFrame for the given cache_data_table function regardless of age.

    Intended as a fallback when a fresh download fails: it returns the last successfully cached
    table even if it's older than the normal 5-day freshness window, or None if no cache exists.
    Which arguments produced that table is deliberately ignored, so a source whose cache key
    includes a version (dbNSFP) can still fall back to the previous version's table on the first
    run after a version bump, rather than failing the whole pipeline.
    The returned table is tagged with FROM_CACHE_ATTR so cache_data_table won't write it back out.
    """
    cache_file_paths = glob.glob(os.path.join(CACHE_DIR, re.sub("^get_", "", function_name) + ".*.tsv.gz"))
    if not cache_file_paths:
        return None
    df = pd.read_table(max(cache_file_paths, key=os.path.getmtime))
    df.attrs[FROM_CACHE_ATTR] = True
    return df


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

        # save result to cache, unless the function fell back to the stale cached table (see
        # read_cached_table). Re-saving that table would bump its mtime and restart the 5-day
        # freshness window above, so it would never expire and the download would never be retried.
        if not df.attrs.get(FROM_CACHE_ATTR):
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