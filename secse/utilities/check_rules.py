import sqlite3
import re
import csv
import os
import sys

def quote_ident(name: str) -> str:
    """用 SQLite 规则安全转义标识符（表名/列名）。"""
    return '"' + name.replace('"', '""') + '"'

def check_smarts(db_file: str, out_file: str):
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        # 写表头
        writer.writerow(["Table", "Rule ID", "SMARTS", "Left Tags", "Right Tags", "Only Left", "Only Right"])

        # 获取所有表和视图
        cur.execute("SELECT name FROM sqlite_master WHERE type IN ('table','view');")
        tables = [r[0] for r in cur.fetchall()]

        pattern = re.compile(r":(\d+)")  # 提取 :数字

        for t in tables:
            qt = quote_ident(t)
            try:
                # 获取表结构
                cur.execute(f"PRAGMA table_info({qt});")
                cols = [row[1] for row in cur.fetchall()]
                lower_cols = [c.lower() for c in cols]

                # 找 smarts 列
                if "smarts" not in lower_cols:
                    continue
                smarts_col = cols[lower_cols.index("smarts")]
                qsmarts = quote_ident(smarts_col)

                # 找 "Rule ID" 列（大小写不敏感）
                id_col = None
                for c in cols:
                    if c.lower().replace(" ", "") in ["ruleid", "rule_id"]:
                        id_col = c
                        break
                if not id_col:
                    continue
                qid = quote_ident(id_col)

                # 查询包含 >> 的行
                cur.execute(f"SELECT {qid}, {qsmarts} FROM {qt} WHERE {qsmarts} LIKE '%>>%';")
                for rid, smarts in cur.fetchall():
                    if not smarts or ">>" not in smarts:
                        continue
                    left, right = smarts.split(">>", 1)

                    left_tags = set(pattern.findall(left))
                    right_tags = set(pattern.findall(right))

                    # 左右标签集合不一致
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
        print("用法: python check_rules.py <input.db> <output.csv>")
        sys.exit(1)

    db_file = sys.argv[1]
    out_file = sys.argv[2]

    if not os.path.exists(db_file):
        print(f"错误: 输入数据库文件不存在 -> {db_file}")
        sys.exit(1)

    try:
        check_smarts(db_file, out_file)
        print(f"✅ 检查完成，结果已保存到 {out_file}")
    except Exception as e:
        print(f"❌ 处理失败: {e}")
        sys.exit(1)

