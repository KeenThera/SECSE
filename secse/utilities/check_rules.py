import sqlite3
import re
import csv
import os
import sys

def quote_ident(name: str) -> str:
    """Quote a table or column identifier using SQLite escaping rules."""
    return '"' + name.replace('"', '""') + '"'

def check_smarts(db_file: str, out_file: str):
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        # Write the header.
        writer.writerow(["Table", "Rule ID", "SMARTS", "Left Tags", "Right Tags", "Only Left", "Only Right"])

        # Get all tables and views.
        cur.execute("SELECT name FROM sqlite_master WHERE type IN ('table','view');")
        tables = [r[0] for r in cur.fetchall()]

        pattern = re.compile(r":(\d+)")  # Extract atom map numbers.

        for t in tables:
            qt = quote_ident(t)
            try:
                # Get the table columns.
                cur.execute(f"PRAGMA table_info({qt});")
                cols = [row[1] for row in cur.fetchall()]
                lower_cols = [c.lower() for c in cols]

                # Find the SMARTS column.
                if "smarts" not in lower_cols:
                    continue
                smarts_col = cols[lower_cols.index("smarts")]
                qsmarts = quote_ident(smarts_col)

                # Find the Rule ID column, ignoring case.
                id_col = None
                for c in cols:
                    if c.lower().replace(" ", "") in ["ruleid", "rule_id"]:
                        id_col = c
                        break
                if not id_col:
                    continue
                qid = quote_ident(id_col)

                # Select reaction SMARTS containing >>.
                cur.execute(f"SELECT {qid}, {qsmarts} FROM {qt} WHERE {qsmarts} LIKE '%>>%';")
                for rid, smarts in cur.fetchall():
                    if not smarts or ">>" not in smarts:
                        continue
                    left, right = smarts.split(">>", 1)

                    left_tags = set(pattern.findall(left))
                    right_tags = set(pattern.findall(right))

                    # Report atom map numbers that differ between sides.
                    if left_tags != right_tags:
                        only_left = sorted(left_tags - right_tags, key=int)
                        only_right = sorted(right_tags - left_tags, key=int)
                        writer.writerow([
                            t,
                            rid,
                            smarts,
                            " ".join(sorted(left_tags, key=int)),
                            " ".join(sorted(right_tags, key=int)),
                            " ".join(only_left),
                            " ".join(only_right)
                        ])

            except Exception:
                continue

    conn.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_rules.py <input.db> <output.csv>")
        sys.exit(1)

    db_file = sys.argv[1]
    out_file = sys.argv[2]

    if not os.path.exists(db_file):
        print(f"Error: input database does not exist: {db_file}")
        sys.exit(1)

    try:
        check_smarts(db_file, out_file)
        print(f"Check complete. Results saved to {out_file}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

