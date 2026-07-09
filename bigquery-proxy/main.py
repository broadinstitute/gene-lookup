import datetime
from flask import jsonify
import gzip
import json
import os
import re
import uuid

from functions_framework import create_app
import google.auth
import google.auth.transport.requests
from google.cloud import bigquery, storage

# Cache the SA email + signing-capable credentials across invocations of the
# warm Cloud Function instance.
_SIGNING_CREDENTIALS = None
_SIGNING_SA_EMAIL = None


def _get_signing_credentials():
    """Return (credentials, service_account_email) for generating V4 signed URLs.

    Cloud Functions Gen2 don't ship a private key, so we sign via the IAM
    signBlob API. The function's runtime SA must have
    roles/iam.serviceAccountTokenCreator on itself.

    The credentials object is cached across invocations on a warm instance, but the underlying
    OAuth access token has a ~1h lifetime, so refresh it whenever it is missing or expired.
    """
    global _SIGNING_CREDENTIALS, _SIGNING_SA_EMAIL
    if _SIGNING_CREDENTIALS is None:
        credentials, _ = google.auth.default()
        # Refresh before reading service_account_email: compute-metadata credentials (Cloud Functions/Run)
        # report service_account_email as the literal "default" until the first refresh populates the real
        # address, so reading it beforehand would cache the wrong signer identity. Treat "default" as
        # unresolved so the K_SERVICE_ACCOUNT / FUNCTION_IDENTITY fallbacks still apply.
        if not credentials.valid:
            credentials.refresh(google.auth.transport.requests.Request())
        sa_email = getattr(credentials, "service_account_email", None)
        if sa_email == "default":
            sa_email = None
        sa_email = sa_email \
            or os.getenv("K_SERVICE_ACCOUNT") \
            or os.getenv("FUNCTION_IDENTITY")
        if not sa_email:
            # user ADC (e.g. `gcloud auth application-default login`) has no service_account_email, so V4
            # signing via IAM signBlob can't work. Fail with an actionable message instead of the opaque
            # AttributeError/ValueError generate_signed_url would otherwise raise. Resolve into locals and
            # raise BEFORE assigning the module globals, so a failed first call leaves the cache empty and a
            # later call (e.g. after the user sets FUNCTION_IDENTITY) re-attempts resolution.
            raise RuntimeError(
                "Cannot generate a V4 signed URL: no service-account identity is available. This path "
                "requires the Cloud Function runtime service account (or the K_SERVICE_ACCOUNT / "
                "FUNCTION_IDENTITY env var). User credentials from `gcloud auth application-default login` "
                "cannot sign; set FUNCTION_IDENTITY to a service account email for local testing."
            )
        _SIGNING_CREDENTIALS = credentials
        _SIGNING_SA_EMAIL = sa_email
    if not _SIGNING_CREDENTIALS.valid:
        _SIGNING_CREDENTIALS.refresh(google.auth.transport.requests.Request())
    return _SIGNING_CREDENTIALS, _SIGNING_SA_EMAIL

BIGQUERY_PROJECT = "cmg-analysis"
BIGQUERY_DATASET = "gene_lookup"
BIGQUERY_TABLE = "combined_gene_disease_association_table"

# Columns that must never be exported to file (OMIM's licensed content). Mirrors allowExport:False in
# website/global_constants.py. Enforced server-side because the client-side export allowlist can be
# bypassed by POSTing raw SQL (e.g. `SELECT *`) directly to this endpoint.
NON_EXPORTABLE_COLUMNS = {
    "OMIM_mim_number",
    "OMIM_phenotype_mim_number",
    "OMIM_inheritance",
    "OMIM_phenotype_description",
}


def _all_schema_field_names(fields):
    """Yield every field name in a BigQuery schema, recursing into RECORD/STRUCT subfields — so a
    non-exportable column nested inside a `SELECT AS STRUCT ...` (whose top-level field is the struct)
    can't hide from the export-column check."""
    for field in fields:
        yield field.name
        if field.fields:
            yield from _all_schema_field_names(field.fields)


def return_query_results(client, result_table, start_index=0, page_size=100):
    """Return query results from BigQuery.

    Args:
        client: The BigQuery client
        result_table: The BigQuery table to return results from
        start_index: The index of the first row to return (default: 0)
        page_size: The number of rows to return (default: 100)

    Returns:
        A tuple of (response_body, response_headers)
    """

    try:
        row_iterator = client.list_rows(result_table, max_results=page_size, start_index=start_index)

        response_dict = {
            "rows": [dict(row) for row in row_iterator],
            "current_page": start_index // page_size + 1,
            "total_pages": (result_table.num_rows + page_size - 1) // page_size if result_table.num_rows > 0 else 0,
            "total_results": result_table.num_rows,
            "page_size": page_size,
            "current_page_start_index": start_index,
            "next_page_start_index": start_index + page_size if start_index + page_size < result_table.num_rows else None,
        }

        return jsonify(response_dict), 200

    except Exception as e:
        print(f"ERROR: Unexpected error: {str(e)}")
        response = jsonify({"error": f"Unexpected error: {e}"})
        return response, 500


def export_to_file(client, result_table, export_to_file_format):
    """Export query results to GCS in the specified format.

    Args:
        client: The BigQuery client
        result_table: The BigQuery table to export results from
        export_to_file_format: The format to export to. Must be one of: 'TSV', 'BED', 'JSON'

    Returns:
        A tuple of (response_body, response_headers)
    """

    if export_to_file_format not in ("TSV", "BED", "JSON"):
        response_dict = {"error": f"Invalid export_to_file_format param: {export_to_file_format}. It must be one of: TSV, BED, JSON"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400

    filename_prefix = uuid.uuid4().hex.upper()
    bucket_name = "gene-lookup-export"
    uri = f'gs://{bucket_name}/{filename_prefix}_*.{export_to_file_format.lower()}.gz'
    if export_to_file_format == "JSON":
        destination_format = bigquery.DestinationFormat.NEWLINE_DELIMITED_JSON
        print_header = None
        field_delimiter = None
    else:
        destination_format = bigquery.DestinationFormat.CSV
        print_header = export_to_file_format == "TSV"
        field_delimiter = '\t'

    try:
        job_config = bigquery.ExtractJobConfig(
            compression=bigquery.Compression.GZIP,
            print_header=print_header,
            field_delimiter=field_delimiter,
            destination_format=destination_format,
        )

        job = client.extract_table(result_table, destination_uris=uri, job_config=job_config)
        job.result()  # Wait for job to complete

        # Get the list of one or more result files
        storage_client = storage.Client()
        bucket = storage_client.bucket(bucket_name)
        blobs = list(bucket.list_blobs(prefix=filename_prefix))
        if len(blobs) == 0:
            response_dict = {"error": f"Export to file did not complete successfully."}
            print(f"ERROR: {response_dict['error']}")
            return jsonify(response_dict), 400

        if export_to_file_format == "JSON" and len(blobs) > 1:
            # Multi-shard exports remain NDJSON; use accurate extension
            export_to_file_format = "NDJSON"

        public_urls = []
        for i, output_blob in enumerate(blobs):
            print(f"DEBUG: Processing output file {i+1} of {len(blobs)}: gs://{bucket_name}/{output_blob.name}")

            if export_to_file_format == "JSON" and len(blobs) == 1:
                temp_file = f"/tmp/{filename_prefix}_original.json.gz"
                new_filename = f"{filename_prefix}_converted.json.gz"
                new_temp_file = f"/tmp/{new_filename}"

                try:
                    # Download the original file
                    output_blob.download_to_filename(temp_file)

                    # Read and convert the newline-delimited JSON to a proper JSON array
                    with gzip.open(temp_file, 'rt') as f, gzip.open(new_temp_file, 'wt') as f_out:
                        f_out.write('[')
                        for j, line in enumerate(f):
                            row = json.loads(line)
                            if j > 0:
                                f_out.write(', ')
                            f_out.write(json.dumps(row, indent=4))
                        f_out.write(']\n')

                    # Upload the converted file to GCS
                    new_blob = bucket.blob(new_filename)
                    new_blob.upload_from_filename(new_temp_file)
                    output_blob.delete()

                    # Update the blob reference to the new file
                    output_blob = new_blob
                finally:
                    # Always remove the local /tmp temp files, even on error, so they don't accumulate
                    # on the warm Cloud Function instance and eventually exhaust its /tmp ramdisk.
                    for f in (temp_file, new_temp_file):
                        try:
                            os.remove(f)
                        except FileNotFoundError:
                            pass

            # Set the Content-Disposition metadata to specify the download filename
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            shard_suffix = f".shard_{i+1:02d}_of_{len(blobs):02d}" if len(blobs) > 1 else ""
            output_filename = f"gene_disease_associations{shard_suffix}.{result_table.num_rows}_genes.{timestamp}.{export_to_file_format.lower()}.gz"
            output_blob.content_type = 'application/gzip'
            output_blob.content_disposition = f'attachment; filename="{output_filename}"'
            output_blob.patch()

            credentials, sa_email = _get_signing_credentials()
            signed_url = output_blob.generate_signed_url(
                version="v4",
                expiration=datetime.timedelta(hours=1),
                method="GET",
                service_account_email=sa_email,
                access_token=credentials.token,
            )
            public_urls.append(signed_url)

        return jsonify({
            'status': 'success',
            'total_results': result_table.num_rows,
            'public_urls': public_urls,
        }), 200

    except Exception as e:
        print(f"ERROR: Unexpected error: {str(e)}")
        return jsonify({"error": f"Unexpected error: {e}"}), 500



def query_gene_lookup_db(request):
    """Cloud Function for forwarding SQL queries to BigQuery.

    Args:
        request (flask.Request): The request object.

        Parameters:
            sql (str): The SQL query to execute.
            start_index (int): (optional)The index of the first row to return (default: 0).
            page_size (int): (optional) The number of rows to return (default: 100).

    Returns:
        flask.Response: The response object.

    """

    response_headers = {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Methods': 'POST',
        'Access-Control-Allow-Headers': 'Content-Type',
    }

    # set CORS headers for the preflight request
    if request.method == 'OPTIONS':
        return '', 204, response_headers

    # validate request body
    try:
        data = request.get_json()
    except Exception as e:
        response_dict = {"error": f"Failed to parse JSON: {request.get_data(as_text=True)}"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    if data is None:
        response_dict = {"error": f"Request body is not valid JSON: {request.get_data(as_text=True)}"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    sql = data.get("sql")
    if not sql:
        response_dict = {"error": f"Missing SQL query in JSON: {request.get_data(as_text=True)}"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    sql = str(sql).strip()

    # Structural validation runs against a normalized copy of the query with (a) single- and double-quoted
    # string literals blanked out — so user search text (which may contain words like "drop", semicolons, or
    # dotted values) can't trip the keyword/table checks below — (b) backtick identifier-quotes removed, so
    # per-identifier quoting like `proj`.`ds`.`tbl` can't smuggle a table reference past the checks, and
    # (c) whitespace around the dot operator collapsed, so `proj . ds . tbl` can't evade the checks.
    # BigQuery accepts both '...' and "..." string literals, so both must be blanked.
    sql_norm = re.sub(r"'(?:[^'\\]|\\.|'')*'", "''", sql)
    sql_norm = re.sub(r'"(?:[^"\\]|\\.|"")*"', '""', sql_norm)
    sql_norm = sql_norm.replace("`", "")
    sql_norm = re.sub(r"\s*\.\s*", ".", sql_norm)

    approved_table = f"{BIGQUERY_PROJECT}.{BIGQUERY_DATASET}.{BIGQUERY_TABLE}"

    # Reject multi-statement scripts. BigQuery's client.query() executes a semicolon-delimited script,
    # so a valid leading SELECT could otherwise be followed by a side-effecting statement (e.g.
    # `EXPORT DATA`, `MERGE`). Allow at most a single optional trailing ';'. Semicolons inside string
    # literals are already blanked out above, so this only sees statement separators.
    if sql_norm.rstrip().rstrip(";").find(";") != -1:
        response_dict = {"error": "Invalid SQL query: only a single SELECT statement is allowed (no ';'-delimited scripts)."}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    # Validate that the query is a SELECT whose first FROM targets exactly the approved table.
    select_pattern = rf"^SELECT\s+.+\s+FROM\s+{re.escape(approved_table)}\b"
    if not re.search(select_pattern, sql_norm, re.IGNORECASE | re.DOTALL):
        response_dict = {"error": f"Invalid SQL query: {sql}. It must be a SELECT from {approved_table}."}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    # After the approved table in the FROM clause, only a standard trailing clause (WHERE/GROUP/HAVING/
    # QUALIFY/WINDOW/ORDER/LIMIT) or the end of the query may follow — never a comma cross-join, JOIN, or
    # table alias introducing a second row source. This blocks row-amplification like
    # `FROM <table>, UNNEST(GENERATE_ARRAY(1, 1e9))` or `FROM <table>, <table>`, whose output cardinality
    # is not bounded by the maximum_bytes_billed (bytes-scanned) cap. UNNEST inside a WHERE subquery is
    # unaffected because it appears after WHERE (the site uses that pattern for rs/transcript/region search).
    from_tail_match = re.search(
        rf"\bFROM\s+{re.escape(approved_table)}\b(?P<tail>.*)$", sql_norm, re.IGNORECASE | re.DOTALL)
    tail = from_tail_match.group("tail") if from_tail_match else ""
    if not re.match(r"\s*($|(WHERE|GROUP|HAVING|QUALIFY|WINDOW|ORDER|LIMIT)\b)", tail, re.IGNORECASE):
        response_dict = {"error": "Invalid SQL query: only a single-table SELECT with standard WHERE/GROUP/ORDER/LIMIT clauses is allowed (no joins, comma-joins, or table functions in the FROM clause)."}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    # Block SQL keywords that could be used for injection or side effects (UNION, subqueries, DML, DDL,
    # BigQuery statements like EXPORT DATA / LOAD DATA / MERGE, scripting, and IAM changes). Checked on
    # the string-blanked copy so a legitimate search for text like "drop attacks" isn't rejected.
    # GENERATE_ARRAY / GENERATE_DATE_ARRAY / GENERATE_TIMESTAMP_ARRAY are row-generation primitives that
    # can synthesize huge row counts without scanning storage (so the bytes-billed cap doesn't bound them);
    # the site never uses them, so they're forbidden anywhere in the query.
    forbidden_pattern = r"\b(UNION|INSERT|UPDATE|DELETE|DROP|CREATE|ALTER|TRUNCATE|MERGE|EXPORT|LOAD|GRANT|REVOKE|ASSERT|BEGIN|DECLARE|EXEC|EXECUTE|CALL|INTO\s+OUTFILE|INFORMATION_SCHEMA|GENERATE_ARRAY|GENERATE_DATE_ARRAY|GENERATE_TIMESTAMP_ARRAY)\b"
    if re.search(forbidden_pattern, sql_norm, re.IGNORECASE):
        response_dict = {"error": f"SQL query contains forbidden keywords: {sql}"}
        print(f"ERROR: {response_dict['error']}")
        return jsonify(response_dict), 400, response_headers

    try:
        client = bigquery.Client()
        # Cap bytes scanned per query to 1 GB (combined_gene_disease_association_table is ~51 MB,
        # so legitimate queries are well under this; the cap bounds cost-drain attacks even if
        # the SQL gate is bypassed).
        max_bytes_billed = 1 * 1024 ** 3

        # Authoritative table allow-list: dry-run the query so BigQuery itself resolves every table
        # reference (including 2-part `dataset.table` names, wildcard tables, JOINs, and subquery FROMs
        # that regex checks miss), then reject unless every referenced table is exactly the approved one.
        # A dry run neither executes the query nor incurs cost. This runs before the real query, so a
        # query that reads a non-approved table never executes.
        dry_run_job = client.query(
            sql, job_config=bigquery.QueryJobConfig(dry_run=True, use_query_cache=False))
        referenced = dry_run_job.referenced_tables
        if not referenced or any(
                (t.project, t.dataset_id, t.table_id) != (BIGQUERY_PROJECT, BIGQUERY_DATASET, BIGQUERY_TABLE)
                for t in referenced):
            response_dict = {"error": f"Invalid SQL query: it may only reference {approved_table}"}
            print(f"ERROR: {response_dict['error']} (referenced: {[t.path for t in referenced]})")
            return jsonify(response_dict), 400, response_headers

        job_config = bigquery.QueryJobConfig(use_query_cache=True, maximum_bytes_billed=max_bytes_billed)
        job = client.query(sql, job_config=job_config)
        job.result()  # wait for the query to complete
        result_table = client.get_table(job.destination)  # get the temporary destination table object
    except Exception as e:
        print(f"ERROR: BigQuery query failed: {str(e)}")
        return jsonify({"error": f"BigQuery query failed: {e}"}), 500, response_headers

    export_to_file_format = data.get("export_to_file_format")
    if export_to_file_format:
        # Enforce the non-exportable (OMIM-licensed) column policy server-side, since a raw POST to this
        # endpoint bypasses the client-side export allowlist. Two complementary checks:
        #  (1) result-schema names, recursing into nested RECORD/STRUCT subfields — catches `SELECT *`,
        #      verbatim `SELECT OMIM_x`, and `SELECT AS STRUCT *` (where OMIM fields are nested);
        #  (2) SELECT-list text (everything before the approved-table FROM) — catches aliases/expressions
        #      like `SELECT OMIM_x AS y` / `CONCAT(OMIM_x, ...)` that (1) misses because the output field
        #      is renamed. Only the SELECT list is checked, so a non-exportable column used purely in a
        #      WHERE filter (allowed — its values never reach the exported file) is not rejected.
        from_match = re.search(rf"\bFROM\s+{re.escape(approved_table)}\b", sql_norm, re.IGNORECASE)
        select_list = sql_norm[:from_match.start()] if from_match else sql_norm
        schema_field_names = set(_all_schema_field_names(result_table.schema))
        disallowed_columns = sorted({
            name for name in schema_field_names if name in NON_EXPORTABLE_COLUMNS
        } | {
            column for column in NON_EXPORTABLE_COLUMNS
            if re.search(rf"\b{re.escape(column)}\b", select_list, re.IGNORECASE)
        })
        if disallowed_columns:
            response_dict = {"error": f"Export is not permitted for these columns: {', '.join(disallowed_columns)}"}
            print(f"ERROR: {response_dict['error']}")
            return jsonify(response_dict), 400, response_headers

        response_json, response_status_code = export_to_file(
            client, result_table, export_to_file_format)

        return response_json, response_status_code, response_headers

    else:
        # Return a page of results
        min_page_size = 1
        default_page_size = 100
        max_page_size = 10_000

        try:
            start_index = int(data.get("start_index", 0))
        except Exception as e:
            response_dict = {"error": f"Invalid start_index param: {e}"}
            print(f"ERROR: {response_dict['error']}")
            return jsonify(response_dict), 400, response_headers

        try:
            page_size = max(min_page_size, min(max_page_size, int(data.get("page_size", default_page_size))))
        except Exception as e:
            response_dict = {"error": f"Invalid page_size param: {e}"}
            print(f"ERROR: {response_dict['error']}")
            return jsonify(response_dict), 400, response_headers

        response_json, response_status_code = return_query_results(
            client, result_table, start_index, page_size)

        return response_json, response_status_code, response_headers

# This is required for Cloud Functions Gen2
if __name__ == "__main__":
    port = int(os.getenv("PORT", "8080"))
    print(f"Starting server on port {port}")
    app = create_app("query_gene_lookup_db")
    app.run(host="0.0.0.0", port=port)
