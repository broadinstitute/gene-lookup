"""Parses JINJA templates and converts them to static HTML pages.

This script processes all *_page_template.html files in the current directory,
renders them using Jinja2, and outputs the final HTML files to the parent directory.
"""

import glob
import json
import jinja2
import os

from global_constants import BIGQUERY_COLUMNS, GROUP_ORDER, CONSTRAINT_THRESHOLDS, get_column_descriptions, get_column_types, get_custom_filter_columns, get_exportable_columns

jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'))

# Get column descriptions to pass to templates
column_descriptions = get_column_descriptions()

# Get column types (STRING/INTEGER/FLOAT) for SQL value-quoting decisions
column_types = get_column_types()

# Get custom filter columns to pass to templates
custom_filter_columns = get_custom_filter_columns()
custom_filter_columns_json = json.dumps(custom_filter_columns)

# Get exportable columns for the export dialog
exportable_columns = get_exportable_columns()
exportable_columns_json = json.dumps(exportable_columns)

with open("data_last_updated_date.json", "r") as f:
    last_updated_dates = json.load(f)

# Per-source version / release-date labels shown next to each source heading on the gene page. A source
# stays absent from this file until it records a successful download, and an absent source shows no label.
# The file itself may not exist yet on a fresh checkout, in which case no source gets a label.
source_versions = json.load(open("source_versions.json")) if os.path.isfile("source_versions.json") else {}

for template_file in glob.glob("*_page_template.html"):
    print(f"Processing {template_file}")

    template = jinja2_env.get_template(os.path.basename(template_file))
    html_content = template.render(
        column_descriptions=column_descriptions,
        column_types=column_types,
        column_groups=GROUP_ORDER,
        custom_filter_columns_json=custom_filter_columns_json,
        exportable_columns_json=exportable_columns_json,
        constraint_thresholds_json=json.dumps(CONSTRAINT_THRESHOLDS),
        data_last_updated_date=last_updated_dates.get("data_last_updated_date", ""),
        source_versions=source_versions,
    )

    output_html_path = os.path.join("../", template_file.replace("_page_template", ""))
    with open(output_html_path, "wt") as f:
        f.write(html_content)

    print(f"Wrote {output_html_path}")
